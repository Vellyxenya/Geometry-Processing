#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/jet.h>
#include <igl/gaussian_curvature.h>
#include <igl/invert_diag.h>
#include <igl/sum.h>
#include <igl/speye.h>
#include <igl/bfs.h>
#include <igl/cotmatrix.h>
#include <igl/principal_curvature.h>
#include <imgui/imgui.h>
/*** insert any libigl headers here ***/
// Custom includes
#include <queue> //queue for custom bfs
#include <set>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #Vx3
Eigen::MatrixXd V;
// Face array, #Fx3
Eigen::MatrixXi F;
//Face normals #Fx3
Eigen::MatrixXd FN;
//Vertex normals #Vx3
Eigen::MatrixXd VN;

// Per-vertex uniform normal array, #Vx3
Eigen::MatrixXd N_uniform;
// Per-vertex area-weighted normal array, #Vx3
Eigen::MatrixXd N_area;
// Per-vertex mean-curvature normal array, #Vx3
Eigen::MatrixXd N_meanCurvature;
// Per-vertex PCA normal array, #Vx3
Eigen::MatrixXd N_PCA;
// Per-vertex quadratic fitted normal array, #Vx3
Eigen::MatrixXd N_quadraticFit;

// Per-vertex mean curvature, #Vx3
Eigen::VectorXd K_mean;
// Per-vertex Gaussian curvature, #Vx3
Eigen::VectorXd K_Gaussian;
// Per-vertex minimal principal curvature, #Vx3
Eigen::VectorXd K_min_principal;
// Per-vertex maximal principal curvature, #Vx3
Eigen::VectorXd K_max_principal;
// Per-vertex color array, #Vx3
Eigen::MatrixXd colors_per_vertex;

// Explicitely smoothed vertex array, #Vx3
Eigen::MatrixXd V_expLap;
// Implicitely smoothed vertex array, #Vx3
Eigen::MatrixXd V_impLap;
// Bilateral smoothed vertex array, #Vx3
Eigen::MatrixXd V_bilateral;

// Custom definitions
int K = 1; //Parametrize the K-Ring
#define _2_PI 6.2831853071795864769
int exp_smoothing_iter = 100;
double lambda = 0.00001;
double delta = 0.001;
double s_c = 0.1;
double s_s = 0.5;
int rho = 1;
int denoise_iter = 5;
bool use_uniform_laplacian = false;
//For curvatures
Eigen::MatrixXd PD_min, PD_max;
Eigen::VectorXd PV_max, PV_min;
Eigen::Vector3d blue(0.2,0.2,0.8);
Eigen::Vector3d red(0.8,0.2,0.2);
//For implicit smoothing
Eigen::SparseMatrix<double> L;

void reorient_normals(Eigen::MatrixXd& normals) {
    Eigen::MatrixXd V_Normals;
    igl::per_vertex_normals(V, F, V_Normals);
    for(int i = 0; i < normals.rows(); i++) {
        double dot = normals.row(i).dot(V_Normals.row(i));
        if(dot < 0) normals.row(i) *= -1;
    }
}

vector<int> KRing(int K, int root, const vector<vector<int>>& adj, Eigen::MatrixXd* neighbors = nullptr) {
    queue<int> q;
    queue<int> depths;
    set<int> explored; //stores the explored nodes
    q.push(root); //enqueue the root
    depths.push(0); //the root is at depth 0
    explored.insert(root); //mark the root as explored
    while(q.size() > 0) {
        int v = q.front();
        q.pop();
        int depth = depths.front();
        depths.pop();
        if(depth < K) {
            for(int j = 0; j < adj[v].size(); j++) { //for each vertex adjacent to i
                if(explored.count(j) <= 0) { //vertex is not explored
                    q.push(adj[v][j]);
                    depths.push(depth + 1); //vertex is 1 layer deeper
                    explored.insert(adj[v][j]); //mark vertex as explored
                }
            }
        }
    }
    vector<int> neighbor_indices(explored.begin(), explored.end());
    if(neighbors != nullptr) {
        neighbors->resize(explored.size(), 3);
        int k = 0;
        for(auto iterator : explored) {
            neighbors->row(k) = V.row(iterator);
            k++;
        }
    }
    return neighbor_indices;
}

void uniform_laplacian(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                       Eigen::SparseMatrix<double>& L) {
    L.resize(V.rows(), V.rows());
    std::vector<Eigen::Triplet<double>> tripletList;
    for(int i = 0; i < F.rows(); i++) {
        for(int j0 = 0; j0 < 3; j0++) {
            //vertex indices
            int j1 = (j0+1) % 3;
            int j2 = (j0+2) % 3;
            int v1 = F(i, j1), v2 = F(i, j2); //vertices
            //For each edge (v1, v2) opposite to current vertex
            //Add a -1 on the corresponding diagonal
            tripletList.push_back(Eigen::Triplet<double>(v1, v1, -1));
            tripletList.push_back(Eigen::Triplet<double>(v2, v2, -1));
            //And add a +1 to indicate this vertex has a neighbor
            tripletList.push_back(Eigen::Triplet<double>(v1, v2, 1));
            tripletList.push_back(Eigen::Triplet<double>(v2, v1, 1));
        }
    }
    L.setFromTriplets(tripletList.begin(), tripletList.end());
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing uniform vertex normals here:
        // store in N_uniform
        N_uniform.setZero(V.rows(), 3);
        igl::per_face_normals(V, F, FN);
        for (int i = 0; i < F.rows(); i++) {
            for (int j = 0; j < 3; j++) {
                N_uniform.row(F(i, j)) += FN.row(i);
            }
        }
        // Reorient the normals
        reorient_normals(N_uniform);
        // Set the viewer normals.
        N_uniform.rowwise().normalize();
        viewer.data().set_normals(N_uniform);
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing area-weighted vertex normals here:
        // store in N_area
        N_area.setZero(V.rows(), 3);
        igl::per_face_normals(V, F, FN);
        Eigen::MatrixXd A;
        igl::doublearea(V, F, A);
        for (int i = 0; i < F.rows(); i++) {
            double area = A(i);
            for (int j = 0; j < 3; j++) {
                N_area.row(F(i, j)) += area * FN.row(i);
            }
        }
        // Reorient the normals
        reorient_normals(N_area);
        // Set the viewer normals.
        N_area.rowwise().normalize();
        viewer.data().set_normals(N_area);
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing mean-curvature vertex normals here:
        // store in N_meanCurvature
        //TODO gives weird 'artifacts' for high-degree vertices (at the center or "hexagons")
        Eigen::SparseMatrix<double> L_c; //Cotangeant Laplacian
        Eigen::SparseMatrix<double> M; //mass matrix
        Eigen::SparseMatrix<double> M_inv; //inverse mass matrix
        igl::cotmatrix(V, F, L_c);
        igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
        igl::invert_diag(M, M_inv);
        Eigen::SparseMatrix<double> L = M_inv * L_c; //Laplacian = area weighted Laplacian
        N_meanCurvature = L * V; //H Matrix
        // Reorient the normals
        reorient_normals(N_meanCurvature);
        // Set the viewer normals.
        N_meanCurvature.rowwise().normalize();
        viewer.data().set_normals(N_meanCurvature);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing PCA vertex normals here:
        // store in N_PCA
        N_PCA.resize(V.rows(), 3);
        vector<vector<int>> A;
        igl::adjacency_list(F, A);
        for (int i = 0; i < V.rows(); i++) {
            Eigen::MatrixXd neighbors;
            KRing(K, i, A, &neighbors);
            Eigen::RowVector3d mu = neighbors.colwise().mean();
            Eigen::MatrixXd X = neighbors.rowwise() - mu;
            Eigen::MatrixXd Cov = X.adjoint() * X;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Cov);
            Eigen::MatrixXd eigenvectors = solver.eigenvectors();
            N_PCA.row(i) = eigenvectors.col(0);
        }
        // Reorient the normals
        reorient_normals(N_PCA);
        // Set the viewer normals.
        N_PCA.rowwise().normalize();
        viewer.data().set_normals(N_PCA);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing quadratic fitted vertex normals here:
        // store in N_quadraticFit
        N_quadraticFit.resize(V.rows(), 3);
        CurvatureCalculator calculator;
        calculator.init(V, F);
        vector<vector<int>> A;
        igl::adjacency_list(F, A);
        for (int i = 0; i < V.rows(); i++) {
            Eigen::MatrixXd neighbors;
            auto neighbor_indices = KRing(K, i, A, &neighbors);
            Eigen::RowVector3d mu = neighbors.colwise().mean();
            Eigen::MatrixXd X = neighbors.rowwise() - mu;
            Eigen::MatrixXd Cov = X.adjoint() * X;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Cov);
            Eigen::MatrixXd eigenvectors = solver.eigenvectors();
            
            //Create the reference frame
            std::vector<Eigen::Vector3d> ref_frame(3);
            ref_frame[0] = eigenvectors.col(1); //u-component
            ref_frame[1] = eigenvectors.col(2); //v-component
            ref_frame[2] = eigenvectors.col(0); //largest eigen vector is z-component

            //Perform quadric fitting
            CurvatureCalculator::Quadric quadric;
            Eigen::Vector3d ref_frame_origin = V.row(i);
            calculator.fitQuadric(ref_frame_origin, ref_frame, neighbor_indices, &quadric);
            
            /**
             * quadric(u, v) = a u^2 + b uv + c v^2 + d u + e v
             * 
             * quadric(u, v) = z => define F(u, v, z) = quadric(u, v) - z  = 0
             * gradF = (du, dv, -1)
             *      du = 2 a u + b v + d
             *      dv = 2 c v + b u + e
             * We evaluate at (u, v) = (0, 0), so the gradient simplifies to:
             * gradF = (d, e, -1)
             * 
             * HOWEVER, the first column of eigenvectors is the z-component,
             * so we get instead: gradF = (-1, d, e)
             **/

            //Normal in local reference frame
            Eigen::RowVector3d local_normal = Eigen::RowVector3d(-1, quadric.d(), quadric.e());
            local_normal.normalize();
            //Transform the normal back to the original reference frame
            N_quadraticFit.row(i) = local_normal * eigenvectors.transpose();
        }
        // Reorient the normals
        reorient_normals(N_quadraticFit);
        // Set the viewer normals.
        N_quadraticFit.rowwise().normalize();
        viewer.data().set_normals(N_quadraticFit);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete mean curvature:
        // store in K_mean
        Eigen::SparseMatrix<double> L_c; //Cotangeant Laplacian
        Eigen::SparseMatrix<double> M; //mass matrix
        Eigen::SparseMatrix<double> M_inv; //inverse mass matrix
        igl::cotmatrix(V, F, L_c);
        igl::massmatrix(V_impLap, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
        igl::invert_diag(M, M_inv);
        Eigen::SparseMatrix<double> L = M_inv * L_c;
        Eigen::MatrixXd H = L * V; //H Matrix
        K_mean = H.rowwise().norm();
        Eigen::MatrixXd V_Normals;
        igl::per_vertex_normals(V, F, V_Normals);
        //curvature sign is given via dot product with the normal
        for(int i = 0; i < H.rows(); i++) {
            double dot = H.row(i).dot(V_Normals.row(i));
            if(dot > 0) K_mean.row(i) *= -1;
        }
        double min = K_mean.minCoeff();
        double max = K_mean.maxCoeff();
        // For visualization, better to normalize the range of K_mean with the maximal and minimal curvatures.
        // store colors in colors_per_vertex
        igl::jet(K_mean, K_mean.minCoeff(), K_mean.maxCoeff(), colors_per_vertex);
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '7') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete Gaussian curvature:
        // store in K_Gaussian
        Eigen::SparseMatrix<double> M;
        Eigen::SparseMatrix<double> M_inv;
        igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
        igl::invert_diag(M, M_inv);
        Eigen::MatrixXd internal_angles;
        igl::internal_angles(V, F, internal_angles);
        K_Gaussian.setConstant(V.rows(), _2_PI);
        for (int i = 0; i < internal_angles.rows(); i++) {
            for (int j = 0; j < 3; j++) {
                K_Gaussian(F(i, j)) -= internal_angles(i, j);
            }
        }
        K_Gaussian = M_inv * K_Gaussian.eval(); //eval() to prevent aliasing
        // For visualization, better to normalize the range of K_Gaussian with the maximal and minimal curvatures.
        // store colors in colors_per_vertex
        double min = K_Gaussian.minCoeff(), max = K_Gaussian.maxCoeff();
        igl::jet(K_Gaussian, min, max, colors_per_vertex);
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '8') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete minimal principal curvature:
        // store in K_min_principal
        K_min_principal.resize(V.rows());
        for (int i = 0; i < V.rows(); i++) {
            double H_i = K_mean(i);
            double G_i = K_Gaussian(i);
            K_min_principal(i) = H_i + sqrt(H_i * H_i - G_i);
        }
        // For visualization, better to normalize the range of K_min_principal with the maximal and minimal curvatures.
        // store colors in colors_per_vertex
        igl::jet(K_min_principal, K_min_principal.minCoeff(), K_min_principal.maxCoeff(), colors_per_vertex);
        
        // Uncomment the code below to draw a blue segment parallel to the minimal curvature direction, 
        const double avg = igl::avg_edge_length(V,F);
        viewer.data().add_edges(V + PD_min*avg, V - PD_min*avg, blue);
        viewer.data().add_edges(V + PD_max*avg, V - PD_max*avg, red);
        
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == '9') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        colors_per_vertex.setZero(V.rows(),3);
        // Add your code for computing the discrete maximal principal curvature:
        // store in K_max_principal
        K_max_principal.resize(V.rows());
        for (int i = 0; i < V.rows(); i++) {
            double H_i = K_mean(i);
            double G_i = K_Gaussian(i);
            K_max_principal(i) = H_i - sqrt(H_i * H_i - G_i);
        }
        // For visualization, better to normalize the range of K_max_principal with the maximal and minimal curvatures
        // store colors in colors_per_vertex
        igl::jet(K_max_principal, K_max_principal.minCoeff(), K_max_principal.maxCoeff(), colors_per_vertex);

        // Uncomment the code below to draw a red segment parallel to the maximal curvature direction
        const double avg = igl::avg_edge_length(V,F);
        viewer.data().add_edges(V + PD_min*avg, V - PD_min*avg, blue);
        viewer.data().add_edges(V + PD_max*avg, V - PD_max*avg, red);
        
        // Set the viewer colors
        viewer.data().set_colors(colors_per_vertex);
    }

    if (key == 'E') {
        // Add your code for computing explicit Laplacian smoothing here:
        // store the smoothed vertices in V_expLap

        Eigen::SparseMatrix<double> I = Eigen::SparseMatrix<double>(V.rows(), V.rows());
        I.setIdentity(); //Sparse Identity matrix for performance
        V_expLap = V;

        /*
        L: area-weighted Laplacian (The Laplacian we actually apply)
        L_c: cotan Laplacian (cotan weights)
        L_u: uniform Laplacian (uniform weights)
        M: mass-matrix that corresponds to the areas

        If we use the cotan Laplacian:
            L = M^(-1) * L_c (depending whether we use uniform or cotan Laplacian)
        If we use the uniform Laplacian:
            L = L_u
        
        The explicit smoothing equation is given by:
        X_(i+1) = (I + t L) X_i
        */
        
        if(use_uniform_laplacian) { //uniform Laplacian
            Eigen::SparseMatrix<double> L_u; //uniform only depends on connectivity
            //so only computed once outside the loop:
            uniform_laplacian(V_expLap, F, L_u);
            //iterate
            for (int i = 0; i < exp_smoothing_iter; i++) {
                V_expLap = (I + lambda * L_u) * V_expLap;
            }
        } else { //cotan Laplacian
            Eigen::SparseMatrix<double> L_c; //cotan Laplacian
            Eigen::SparseMatrix<double> M; //mass matrix
            Eigen::SparseMatrix<double> M_inv; //inverse mass matrix
            for (int i = 0; i < exp_smoothing_iter; i++) {
                igl::cotmatrix(V_expLap, F, L_c);
                igl::massmatrix(V_expLap, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
                igl::invert_diag(M, M_inv);
                L = M_inv * L_c;
                V_expLap = (I + lambda * L) * V_expLap;
            }
        }

        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_expLap, F);
    }

    if (key == 'D'){
        // Implicit smoothing for comparison
        // store the smoothed vertices in V_impLap
        // Taken from https://github.com/libigl/libigl/blob/main/tutorial/205_Laplacian/main.cpp

        // Recompute just mass matrix on each step
        Eigen::SparseMatrix<double> M;
        igl::massmatrix(V_impLap,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
        // Solve (M-delta*L) U = M*U
        const auto & S = (M - delta*L);
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
        assert(solver.info() == Eigen::Success);
        V_impLap = solver.solve(M*V_impLap).eval();
        // Compute centroid and subtract (also important for numerics)
        Eigen::VectorXd dblA;
        igl::doublearea(V,F,dblA);
        double area = 0.5*dblA.sum();
        Eigen::MatrixXd BC;
        igl::barycenter(V_impLap,F,BC);
        Eigen::RowVector3d centroid(0,0,0);
        for(int i = 0;i<BC.rows();i++)
        {
          centroid += 0.5*dblA(i)/area*BC.row(i);
        }
        V_impLap.rowwise() -= centroid;
        // Normalize to unit surface area (important for numerics)
        //V_impLap.array() /= sqrt(area); //Commented out otherwise shrinks and is annoying
        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_impLap, F);
    }

    if (key == 'B') {
        // Add your code for computing bilateral smoothing here:
        // store the smoothed vertices in V_bilateral
        V_bilateral = V;
        vector<vector<int>> A; //adjacency list
        igl::adjacency_list(F, A);
        for (int k = 0; k < denoise_iter; k++) {
            Eigen::MatrixXd V_Normals;
            igl::per_vertex_normals(V_bilateral, F, V_Normals); //Recompute the new normals
            Eigen::MatrixXd V_bilateral_temp = V_bilateral;
            for (int i = 0; i < V.rows(); i++) {
                std::vector<int> neighbors = KRing(rho, i, A);
                double sum = 0;
                double normalizer = 0;
                for (int j = 0; j < neighbors.size(); j++) {
                    auto v_minus_q_i = V_bilateral.row(i) - V_bilateral.row(neighbors[j]);
                    double t = (v_minus_q_i).norm();
                    double h = V_Normals.row(i).dot(v_minus_q_i);
                    double w_c = exp(-t * t / (2 * s_c * s_c));
                    double w_s = exp(-h * h / (2 * s_s * s_s));
                    double wc_ws = w_c * w_s;
                    sum += wc_ws * h;
                    normalizer += wc_ws;
                }
                V_bilateral_temp.row(i) = V_bilateral.row(i) - V_Normals.row(i) * sum / normalizer;
            }
            V_bilateral = V_bilateral_temp;
        }

        // Set the smoothed mesh
        viewer.data().clear();
        viewer.data().set_mesh(V_bilateral, F);
    }

    if (key == 'R') {
        //Reset the implicitly smoothed vertices to their original values
        V_impLap = V;
        viewer.data().clear();
        viewer.data().set_mesh(V_impLap, F);
    }

    return true;
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    if (filename.substr(filename.length() - 4) == ".off")
    {
        igl::readOFF(filename, V, F);
    }
    else if (filename.substr(filename.length() - 4) == ".obj")
    {
        igl::readOBJ(filename, V, F);
    }
    else
    {
        std::cerr << "Extension unknown (must be '.off' or '.obj')\n";
        return false;
    }
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().compute_normals();
    viewer.core().align_camera_center(V, F);
    return true;
}

int main(int argc, char *argv[]) {
    // Show the mesh
    Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    std::string filename;
    if (argc == 2) {
        filename = std::string(argv[1]);
    }
    else {
        filename = std::string("../data/bumpy-cube.obj");
    }
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;

    menu.callback_draw_viewer_menu = [&]() {
        // Draw parent menu content
        menu.draw_viewer_menu();
        // Add new group
        if (ImGui::CollapsingHeader("Normals", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputInt("K-Ring", &K);
        }
        if (ImGui::CollapsingHeader("Explicit Laplacian smoothing", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputInt("Iterations", &exp_smoothing_iter);
            ImGui::InputDouble("lambda", &lambda, 0.0, 0.0, "%.8f");
            ImGui::Checkbox("Uniform Laplacian", &use_uniform_laplacian);
        }
        if (ImGui::CollapsingHeader("Implicit Laplacian smoothing", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputDouble("delta", &delta, 0.0, 0.0, "%.8f");
        }
        if (ImGui::CollapsingHeader("Bilaterial denoising", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::InputDouble("sigma_c", &s_c, 0.0, 0.0, "%.6f");
            ImGui::InputDouble("sigma_s", &s_s, 0.0, 0.0, "%.6f");
            ImGui::InputInt("K-ring", &rho);
            ImGui::InputInt("Denoising iter", &denoise_iter);
        }
    };

    viewer.plugins.push_back(&menu);

    //Initializations
    // Compute Laplace-Beltrami operator: #V by #V (Implicit smoothing)
    igl::cotmatrix(V,F,L);
    V_impLap = V;
    // Compute principal curvature directions
    
    igl::principal_curvature(V, F, PD_max, PD_min, PV_max, PV_min);
    
    viewer.launch();

}
