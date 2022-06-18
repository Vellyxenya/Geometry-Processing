#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>

/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/dijkstra.h>
#include <random>


using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

// Custom initializations
// Random number generation
std::random_device rd;
std::mt19937 rng(rd());
std::uniform_int_distribution<int> uni;
// Strategy to choose the 2 fixed points
enum FixedPointsStrategy {
	BF_GEODESIC = 0, //Brute force geodesic distance using dijkstra
	BF_EUCLIDEAN, //Brute force euclidean distance
	MC_GEODESIC, //Works but ALTERNATING_MC_ITERATIVE_GEODESIC is consistently better
	ITERATIVE_GEODESIC, //The best for ARAP as we want to keep the same fixed points between iterations (no randomness involved here)
	ALTERNATING_MC_ITERATIVE_GEODESIC, //The best for non-ARAP :D
};
FixedPointsStrategy strategy = ALTERNATING_MC_ITERATIVE_GEODESIC;
//distortions
enum Distorsion { None = 0, Conformal, Isometric, Authalic, Dirichlet };
Distorsion distortion = None;
Eigen::MatrixXd colors; 

void Redraw()
{
	viewer.data().clear();

	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(false);

    if(UV.size() != 0)
    {
      viewer.data().set_uv(TextureResolution*UV);
      viewer.data().show_texture = true;
    }
	}
	else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);
	}
	//Set the colors for distortion visualization
	viewer.data().set_colors(colors);
}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}

static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}

void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double> &C, VectorXd &d)
{
	// Convert the list of fixed indices and their fixed positions to a linear system
	// Hint: The matrix C should contain only one non-zero element per row and d should contain the positions in the correct order.
	int P = positions.rows();
	C.resize(P * 2, V.rows() * 2);
	d.resize(P * 2); //Flatten x and y dimensions
	vector<Eigen::Triplet<double>> data;
	for (int i = 0; i < P; i++) {
		d(i) = positions(i, 0);
		d(i + P) = positions(i, 1);
		data.push_back(Eigen::Triplet<double>(i, indices(i), 1));
		data.push_back(Eigen::Triplet<double>(i + P, indices(i) + V.rows(), 1));
	}
	C.setFromTriplets(data.begin(), data.end());
}

void localStepARAP(Eigen::MatrixXd& R,
			   	   const Eigen::SparseMatrix<double>& Dx,
			   	   const Eigen::SparseMatrix<double>& Dy) {
	Eigen::VectorXd Dxu(Dx * UV.col(0)), Dyu(Dy * UV.col(0)),
		Dxv(Dx * UV.col(1)), Dyv(Dy * UV.col(1));
	for (int i = 0; i < F.rows(); i++) { //Iterate over triangles
		Eigen::Matrix2d J, U, S, V;
		J << Dxu[i], Dyu[i],
			 Dxv[i], Dyv[i];
		SSVD2x2(J, U, S, V);
		Eigen::Matrix2d R_J = U * V.transpose();
		R.row(i) << R_J(0, 0), R_J(0, 1), R_J(1, 0), R_J(1, 1);
	}
}

double measure_distortion(const Eigen::Matrix2d& J) {
	double measure = 0;
	Eigen::Matrix2d U, S, V, R;
	switch (distortion) {
	case Conformal:
		measure = (J + J.transpose() - J.trace() * Eigen::Matrix2d::Identity()).norm();
		break;
	case Isometric:
		SSVD2x2(J, U, S, V);
		R = U * V.transpose();
		measure = (J - R).norm();
		break;
	case Authalic:
		measure = pow(J.determinant() - 1, 2);
		break;
	case Dirichlet:
		measure = J.norm();
		break;
	default:
		cout << "No distortion measure provided" << endl;
		break;
	}
	return measure;
}

void visualize_distortion() {
	colors = Eigen::MatrixXd(F.rows(), 3);
	if (distortion == None) {
		colors.rowwise() = Eigen::RowVector3d(1., 1., 1.);
	}
	Eigen::SparseMatrix<double> Dx, Dy;
	computeSurfaceGradientMatrix(Dx, Dy);
	Eigen::VectorXd Dxu(Dx * UV.col(0)), Dyu(Dy * UV.col(0)), Dxv(Dx * UV.col(1)), Dyv(Dy * UV.col(1));
	vector<double> distortions(F.rows());
	for (int i = 0; i < F.rows(); i++) { //Iterate over triangles
		Eigen::Matrix2d J, U, S, V;
		J << Dxu[i], Dyu[i],
			Dxv[i], Dyv[i];
		distortions[i] = measure_distortion(J);
	}
	double max_distortion = *max_element(distortions.begin(), distortions.end());
	for (int i = 0; i < distortions.size(); i++) {
		double intensity = 1 - distortions[i] / max_distortion;
		colors.row(i) = Eigen::RowVector3d(1., intensity, intensity);
	}
}

void computeParameterization(int type)
{
	VectorXi fixed_UV_indices;
	MatrixXd fixed_UV_positions;

	SparseMatrix<double> A;
	VectorXd b;
	Eigen::SparseMatrix<double> C;
	VectorXd d;
	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	if (!freeBoundary) {
		// The boundary vertices should be fixed to positions on the unit disc. Find these position and
		// save them in the #V x 2 matrix fixed_UV_position.
		igl::boundary_loop(F, fixed_UV_indices);
		igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
	} else {
		// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.
		vector<vector<int>> Adj;
		Eigen::VectorXd min_dist;
		igl::adjacency_list(F, Adj);
		int max_dist = 0;
		int arg_i = 0, arg_j = 0;
		int max_idx_i;
		if(strategy == BF_GEODESIC) { //Runs in O(V * E logV) > O(V²)
			//Brute force geodesic distance metric
			//NOTE: This is perfect for small meshes but too slow for large ones
			for (int i = 0; i < V.rows(); i++) {
				Eigen::VectorXi prev;
				igl::dijkstra(i, set<int>(), Adj, min_dist, prev);
				int max_dist_i = min_dist.maxCoeff(&max_idx_i);
				if (max_dist_i > max_dist) {
					max_dist = max_dist_i;
					arg_i = i;
					arg_j = max_idx_i;
					cout << max_dist << " " << arg_i << " " << arg_j << endl;
				}
			}
		} else if(strategy == BF_EUCLIDEAN) { //Runs in O(V²)
			//Brute force euclidean distance metric
			for (int i = 0; i < V.rows(); i++) { 
				for (int j = i + 1; j < V.rows(); j++) { 
					float dist = (V.row(i) - V.row(j)).norm();
					if(dist > max_dist) {
						arg_i = i;
						arg_j = j;
						max_dist = dist;
					}
				}
			}
		} else if(strategy == ITERATIVE_GEODESIC) {
			int start = 0;
			int max_iter = 100;
			Eigen::VectorXi prev;
			for(int i = 0; i < max_iter; i++) {
				igl::dijkstra(start, set<int>(), Adj, min_dist, prev);
				int max_dist_i = min_dist.maxCoeff(&max_idx_i);
				if (max_dist_i > max_dist) {
					max_dist = max_dist_i;
					arg_i = start;
					arg_j = max_idx_i;
					cout << max_dist << " " << arg_i << " " << arg_j << endl;
				}
				//New start is point with max distance to previous start
				//Guaranteed to monotonically increase the distance
				//Convergence rate?
				//Guaranteed to find optimal solution?
				start = max_idx_i;
			}
		} else if(strategy == MC_GEODESIC) {
			int start;
			int max_iter = 100;
			Eigen::VectorXi prev;
			for(int i = 0; i < max_iter; i++) {
				start = uni(rng);
				igl::dijkstra(start, set<int>(), Adj, min_dist, prev);
				int max_dist_i = min_dist.maxCoeff(&max_idx_i);
				if (max_dist_i > max_dist) {
					max_dist = max_dist_i;
					arg_i = start;
					arg_j = max_idx_i;
					cout << max_dist << " " << arg_i << " " << arg_j << endl;
				}
			}
		} else if(strategy == ALTERNATING_MC_ITERATIVE_GEODESIC) {
			int start;
			int max_iter = 100;
			Eigen::VectorXi prev;
			for(int i = 0; i < max_iter; i++) {
				start = i%2 == 0 ? uni(rng) : start; //Alternate between a fresh start and 
													 //Greed based on previous best result
				igl::dijkstra(start, set<int>(), Adj, min_dist, prev);
				int max_dist_i = min_dist.maxCoeff(&max_idx_i);
				if (max_dist_i > max_dist) {
					max_dist = max_dist_i;
					arg_i = start;
					arg_j = max_idx_i;
					cout << max_dist << " " << arg_i << " " << arg_j << endl;
				}
				start = max_idx_i;
			}
		}
		fixed_UV_indices.resize(2);
		fixed_UV_indices << arg_i, arg_j;
		fixed_UV_positions.resize(2, 2);
		fixed_UV_positions << 1, 0, -1, 0;
		cout << "Finished computing the 2 fixed vertices" << endl;
	}

	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);

	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh
		Eigen::SparseMatrix<double> L;
		igl::adjacency_matrix(F, L);
		Eigen::VectorXd sum_vec;
		igl::sum(L, 2, sum_vec); //rowwise sum, dim index seems to start at 1 for some reason
		//Very important to cast for some reason.. Otherwise very obscure crash
		L -= (Eigen::SparseMatrix<double>)sum_vec.asDiagonal(); //Uniform Laplacian

		Eigen::SparseMatrix<double> Zero(V.rows(), V.rows());
		Eigen::SparseMatrix<double> Top = igl::cat(2, L, Zero);
		Eigen::SparseMatrix<double> Bottom = igl::cat(2, Zero, L);
		A = igl::cat(1, Top, Bottom);
		b.setZero(2 * V.rows());
	}

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~
		Eigen::SparseMatrix<double> L_c;
		Eigen::SparseMatrix<double> M;
		Eigen::SparseMatrix<double> M_inv;
		igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		igl::invert_diag(M, M_inv);
		igl::cotmatrix(V, F, L_c);
		//Gives visually strictly similar results without using the mass matrix.
		//Is it needed then?
		Eigen::SparseMatrix<double> L = M_inv * L_c;
		Eigen::SparseMatrix<double> Zero(V.rows(), V.rows());
		Eigen::SparseMatrix<double> Top = igl::cat(2, L, Zero);
		Eigen::SparseMatrix<double> Bottom = igl::cat(2, Zero, L);
		A = igl::cat(1, Top, Bottom);
		b.setZero(2 * V.rows());
	}

	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!
		Eigen::SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(Dx, Dy);
		Eigen::VectorXd double_area;
		igl::doublearea(V, F, double_area);
		//The A_cal = A^T A includes '2' factors everywhere, so no need to divide 'double_area' by 2.
		Eigen::SparseMatrix<double> A_f = (Eigen::SparseMatrix<double>)double_area.asDiagonal();
		Eigen::SparseMatrix<double> DxT(Dx.transpose()), DyT(Dy.transpose());
		Eigen::SparseMatrix<double> TopLeft     = (DxT * A_f * Dx + DyT * A_f * Dy);
		Eigen::SparseMatrix<double> TopRight    = (DyT * A_f * Dx - DxT * A_f * Dy);
		Eigen::SparseMatrix<double> BottomLeft  = (DxT * A_f * Dy - DyT * A_f * Dx);
		Eigen::SparseMatrix<double> BottomRight = (DxT * A_f * Dx + DyT * A_f * Dy);
		Eigen::SparseMatrix<double> Top    = igl::cat(2, TopLeft, TopRight);
		Eigen::SparseMatrix<double> Bottom = igl::cat(2, BottomLeft, BottomRight);
		A = igl::cat(1, Top, Bottom);
		b.setZero(2 * V.rows());
		//Sanity Check: Indeed with fixed boundary, both '2' and '3' give same parametrization
	}

	if (type == '4') {
		//NOTE 1: Compute LSCM parametrization first by pressing '3'
		//NOTE 2: Press '4' repeatedly until it converges
		
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices
	
		Eigen::SparseMatrix<double> Dx, Dy;
		computeSurfaceGradientMatrix(Dx, Dy);
		Eigen::MatrixXd R(F.rows(), 4);

		//Find closest rotation matrices
		localStepARAP(R, Dx, Dy);

		Eigen::VectorXd double_area;
		igl::doublearea(V, F, double_area);
		Eigen::SparseMatrix<double> A_f = (Eigen::SparseMatrix<double>)(double_area/2).asDiagonal();
		Eigen::SparseMatrix<double> DxT(Dx.transpose()), DyT(Dy.transpose());
		Eigen::SparseMatrix<double> Zero(V.rows(), V.rows());
		Eigen::SparseMatrix<double> L = (DxT * A_f * Dx + DyT * A_f * Dy);
		Eigen::SparseMatrix<double> Top = igl::cat(2, L, Zero);
		Eigen::SparseMatrix<double> Bottom = igl::cat(2, Zero, L);
		A = igl::cat(1, Top, Bottom);
		Eigen::VectorXd b1 = DxT * A_f * R.col(0) + DyT * A_f * R.col(1); //R.col(0) = R11, R.col(1) = R12
		Eigen::VectorXd b2 = DxT * A_f * R.col(2) + DyT * A_f * R.col(3); //R.col(2) = R21, R.col(3) = R22
		b = igl::cat(1, b1, b2);
	}

	// Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet
	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail
	Eigen::SparseLU<Eigen::SparseMatrix<double>> sparse_LU_solver;
	/* E x = e */
	//Setup Matrix and vector
	Eigen::SparseMatrix<double> E;
	Eigen::SparseMatrix<double> Zero(d.size(), d.size());
	Eigen::SparseMatrix<double> CT = C.transpose();
	Eigen::SparseMatrix<double> Top = igl::cat(2, A, CT);
	Eigen::SparseMatrix<double> Bottom = igl::cat(2, C, Zero);
	E = igl::cat(1, Top, Bottom);
	Eigen::VectorXd e;
	e.resize(b.size() + d.size());
	e << b, d;
	//Solve
	sparse_LU_solver.analyzePattern(E);
	sparse_LU_solver.factorize(E);
	Eigen::VectorXd u_v_lambda = sparse_LU_solver.solve(e); //lambda will be discarded

	// The solver will output a vector
	UV.resize(V.rows(), 2);
	UV.col(0) = u_v_lambda.block(0, 0, V.rows(), 1); //u
	UV.col(1) = u_v_lambda.block(V.rows(), 0, V.rows(), 1); //v
}

bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
	case '4':
		computeParameterization(key);
		break;
	case '5':
		visualize_distortion();
		break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
    if(showingUV)
    {
      temp2D = viewer.core();
      viewer.core() = temp3D;
      showingUV = false;
    }
    else
    {
      if(UV.rows() > 0)
      {
        temp3D = viewer.core();
        viewer.core() = temp2D;
        showingUV = true;
      }
      else { std::cout << "ERROR ! No valid parameterization\n"; }
    }
    break;
	}
	Redraw();
	return true;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  Redraw();
  viewer.core().align_camera_center(V);
  showingUV = false;

  return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core();
	temp2D = viewer.core();
	temp2D.orthographic = true;

	return false;
}

int main(int argc,char *argv[]) {
	if(argc != 2) {
		cout << "Usage ex4_bin <mesh.off/obj>" << endl;
		load_mesh("../data/cathead.obj");
	}
	else
	{
		// Read points and normals
		load_mesh(argv[1]);
	}

	//Initialize uniform sampler
	uni = std::uniform_int_distribution<int>(0, V.rows()-1);

	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parmatrization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);
			const char* names[] = {"Brute Force Geodesic", "Brute Force Euclidean",
				"Monte Carlo Geodesic", "Iterative Geodesic", "Alternating MC-Iterative Geodesic"
			};
			ImGui::Combo("Fixed indices strategy", (int*)&strategy, names, 5);

			const char* distortion_names[] = {"None", "Conformal",
				"Isometric", "Authalic", "Dirichlet"
			};
			ImGui::Combo("Distortion", (int*)&distortion, distortion_names, 5);
		}
	};

	viewer.callback_key_pressed = callback_key_pressed;
	viewer.callback_mouse_move = callback_mouse_move;
	viewer.callback_init = callback_init;

	viewer.launch();
}
