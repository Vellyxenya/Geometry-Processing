#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/vertex_components.h>
#include <chrono>

using namespace std::chrono;
using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Normals evaluated via PCA method, #P x3
Eigen::MatrixXd NP;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Parameter: degree of the polynomial
int polyDegree = 0;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
double wendlandRadius = 0.1;

// Parameter: grid resolution
//int resolution = 20;

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

/*
Own definitions
*/
#pragma region
// Predefined colors:
const Eigen::RowVector3d RED  (1, 0, 0);
const Eigen::RowVector3d GREEN(0, 1, 0);
const Eigen::RowVector3d BLUE (0, 0, 1);
// Some global variables
Eigen::RowVector3d offset; // origin of the 3d grid
double EPSILON; // original epsilon
double cell_size; // cell side length
int rX, rY, rZ; // resolutions in each axis
vector<vector<vector<vector<int>>>> spatial_index; //data structure to accelerate neighbor queries
double inflation_ratio = 1.1; //Used so that the grid bounds are not too close to the mesh
bool pca_bb = false;
string file_name; // input file
double pca_radius = 10;
int resolutionX = 20;
int resolutionY = 20;
int resolutionZ = 20;
bool connected_components = false;
#pragma endregion

/*
Own functions
*/
#pragma region 
void create_spatial_index() {
    Eigen::RowVector3d min = P.colwise().minCoeff(), max = P.colwise().maxCoeff();
    Eigen::RowVector3d diff = max - min;
    cell_size = 2 * EPSILON;
    rX = diff[0] / cell_size + 1; // crashes without the +1
    rY = diff[1] / cell_size + 1;
    rZ = diff[2] / cell_size + 1;
    offset = min;
    //cout << "cell_size: " << cell_size << endl;
    //cout << "dims: " << rX << " " << rY << " " << rZ << endl;
    spatial_index = vector<vector<vector<vector<int>>>>(rX, vector<vector<vector<int>>>(rY, vector<vector<int>>(rZ, vector<int>())));
    for (int i = 0; i < P.rows(); i++) {
        Eigen::RowVector3d p = (P.row(i) - offset) / cell_size;
        int indexX = (int)p.x();
        int indexY = (int)p.y();
        int indexZ = (int)p.z();
        /*if (indexX >= rX || indexY >= rY || indexZ >= rZ || indexX < 0 || indexY < 0 || indexZ < 0) {
            cout << indexX << " " << indexY << " " << indexZ << endl;
            cout << p << endl;
        }*/
        spatial_index[indexX][indexY][indexZ].push_back(i);
    }
}

bool brute_force_proximity(const Eigen::RowVector3d& point, int index) {
    double squared_norm = (point - P.row(index)).squaredNorm();
    for (int i = 0; i < P.rows(); i++) {
        double squared_norm_i = (point - P.row(i)).squaredNorm();
        if (squared_norm_i < squared_norm) return false;
    }
    return true;
}

bool spatial_index_proximity(const Eigen::RowVector3d& point, int index) {
    double squared_norm = (point - P.row(index)).squaredNorm();
    Eigen::RowVectorXd p = (point - offset) / cell_size;
    int px = (int)p.x(), py = (int)p.y(), pz = (int)p.z();
    for (int i = max(px - 1, 0); i <= min(px + 1, rX - 1); i++) {
        for (int j = max(py - 1, 0); j <= min(py + 1, rY - 1); j++) {
            for (int k = max(pz - 1, 0); k <= min(pz + 1, rZ - 1); k++) {
                for (int p_i : spatial_index[i][j][k]) {
                    double squared_norm_i = (point - P.row(p_i)).squaredNorm();
                    if (squared_norm_i < squared_norm) return false;
                }
            }
        }
    }
    return true;
}

bool original_is_nearest(const Eigen::RowVector3d& point, int index) {
    //return brute_force_proximity(point, index); //For performance measurement
    return spatial_index_proximity(point, index);
}

int get_basis_dim() {
    if (polyDegree > 2) {
        cout << "Only support polyDegree up to 2, setting polyDegree to 2" << endl;
        polyDegree = 2;
    }
    int dims[] = { 1, 4, 10 };
    return dims[polyDegree];
}

Eigen::RowVectorXd get_basis_vector(Eigen::RowVector3d point, int dim) {
    Eigen::RowVectorXd basis;
    basis.resize(10); //TODO necessary?
    double x = point.x(), y = point.y(), z = point.z();
    basis << 1, x, y, z, x * x, y * y, z * z, x * y, x * z, y * z;
    return basis.block(0, 0, 1, dim);
}

double wendland(double r) {
    return pow(1 - r/wendlandRadius, 4) * (4 * r / wendlandRadius + 1);
}

double S_k(Eigen::RowVectorXd x, int i) {
    double normal_norm = N.row(i).norm();
    Eigen::RowVectorXd normalized_normal = N.row(i) / normal_norm;
    return constrained_values(i) + (x - P.row(i)).dot(normalized_normal);
}

vector<int> get_neighbors(const Eigen::RowVector3d& point, bool only_originals = false, float search_radius = -1) {
    Eigen::RowVectorXd p = (point - offset) / cell_size;
    vector<int> neighbors;
    search_radius = (search_radius == -1)? wendlandRadius : search_radius;
    double radiusSquared = search_radius * search_radius;
    int delta = ceil(search_radius / cell_size); // +1;
    //cout << "delta: " << delta << endl;
    int px = (int)p.x(), py = (int)p.y(), pz = (int)p.z();
    for (int i = max(px - delta, 0); i <= min(px + delta, rX - 1); i++) {
        for (int j = max(py - delta, 0); j <= min(py + delta, rY - 1); j++) {
            for (int k = max(pz - delta, 0); k <= min(pz + delta, rZ - 1); k++) {
                for (int p_i_0 : spatial_index[i][j][k]) {
                    double squared_norm_i = (point - constrained_points.row(p_i_0)).squaredNorm();
                    if (squared_norm_i < radiusSquared) {
                        neighbors.push_back(p_i_0);
                        if(!only_originals) {
                            neighbors.push_back(p_i_0 + (int)P.rows());
                            neighbors.push_back(p_i_0 + 2 * (int)P.rows());
                        }
                    }
                    /*vector<int> ps = {p_i_0, p_i_0 + (int)P.rows(), p_i_0 + 2 * (int)P.rows()};
                    for (int p_i : ps) {
                        double squared_norm_i = (point - constrained_points.row(p_i)).squaredNorm();
                        if (squared_norm_i < wendlandRadiusSquared)
                            neighbors.push_back(p_i);
                    }*/
                }
            }
        }
    }
    /*for (int i = 0; i < constrained_points.rows(); i++) {
        double squared_norm_i = (point - constrained_points.row(i)).squaredNorm();
        if (squared_norm_i < wendlandRadiusSquared)
            neighbors.push_back(i);
    }*/

    return neighbors;
}

vector<int> get_neighbors_brute(const Eigen::RowVector3d& point) {
    double squared_radius = wendlandRadius * wendlandRadius;
    vector<int> neighbors;
    for(int i = 0; i < P.rows(); i++) {
        double squared_norm = (P.row(i) - point).squaredNorm();
        if(squared_norm <= squared_radius) {
            neighbors.push_back(i);
            neighbors.push_back(i + (int)P.rows());
            neighbors.push_back(i + 2 * (int)P.rows());
        }
    }
    return neighbors;
}

void run_connected_components(Eigen::MatrixXd& V_, Eigen::MatrixXi& F_) {
    Eigen::VectorXi cid;
    igl::vertex_components(F, cid);
    int nb_components = cid.maxCoeff() + 1;
    std::vector<int> nb_vertices_per_component(nb_components, 0);
    for (int i = 0; i < cid.size(); i++) {
        nb_vertices_per_component[cid(i)]++;
    }
    vector<int>::iterator it = max_element(nb_vertices_per_component.begin(), nb_vertices_per_component.end());
    int largest_component_id = distance(nb_vertices_per_component.begin(), it);
    int k = 0;
    std::map<int, int> map;
    Eigen::MatrixXd V_new;
    Eigen::MatrixXi F_new;
    V_new.resize(V.rows(), 3);
    F_new.resize(F.rows(), 3);
    for (int i = 0; i < cid.size(); i++) {
        if (cid(i) == largest_component_id) {
            V_new.row(k) = V.row(i);
            map[i] = k;
            k++;
        }
    }
    V_ = V_new.block(0, 0, k, 3);
    k = 0;
    for (int i = 0; i < F.rows(); i++) {
        if (map.count(F(i, 0)) > 0) { //If first vertex is in the new vertex set (meaning the face must be kept)
            F_new.row(k) = Eigen::RowVector3i(map[F(i, 0)], map[F(i, 1)], map[F(i, 2)]);
            k++;
        }
    }
    F_ = F_new.block(0, 0, k, 3);
}
#pragma endregion


// Functions
void createGrid();
void evaluateImplicitFunc();
void evaluateImplicitFunc_PolygonSoup();
void getLines();
void pcaNormal();
bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers);

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid()
{
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines.resize(0, 6);
    grid_values.resize(0);
    V.resize(0, 3);
    F.resize(0, 3);
    FN.resize(0, 3);

    Eigen::MatrixXd W = Eigen::MatrixXd::Identity(3, 3); //If no PCA, the transformation is the identity
    Eigen::RowVector3d mu = P.rowwise().mean();
    Eigen::MatrixXd X = P.rowwise() - mu;

    if (pca_bb) {
        //Perform PCA on the data points
        Eigen::MatrixXd Cov = X.adjoint() * X;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Cov);
        Eigen::MatrixXd eigenvectors = solver.eigenvectors();
        W = eigenvectors.determinant() < 0 ? -eigenvectors : eigenvectors;
    }
    
    //Transform the data
    Eigen::MatrixXd T = X * W;

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = T.colwise().minCoeff();
    bb_max = T.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d center = (bb_max + bb_min) / 2;
    bb_min = center + (bb_min - center) * inflation_ratio;
    bb_max = center + (bb_max - center) * inflation_ratio;
    Eigen::RowVector3d dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolutionX - 1);
    const double dy = dim[1] / (double)(resolutionY - 1);
    const double dz = dim[2] / (double)(resolutionZ - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolutionX * resolutionY * resolutionZ, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolutionX; ++x)
    {
        for (unsigned int y = 0; y < resolutionY; ++y)
        {
            for (unsigned int z = 0; z < resolutionZ; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolutionX * (y + resolutionY * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }

    //Transform the grid_points back to original coordinate system
    grid_points = grid_points * W.inverse();
    grid_points = grid_points.rowwise() + mu;
}

// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc()
{
    grid_values.resize(resolutionX * resolutionY * resolutionZ);
    Eigen::MatrixXd B;
    Eigen::VectorXd W;
    Eigen::VectorXd F;
    int basis_dim = get_basis_dim();
    auto start = high_resolution_clock::now();
    for (int i = 0; i < grid_points.rows(); i++) {
        //Notations follow those from the course
        const auto x = grid_points.row(i);
        vector<int> neighbors = get_neighbors(x);
        //vector<int> neighbors = get_neighbors_brute(x); //For performance measurement
        if (neighbors.empty()) grid_values(i) = INT_MAX; //Assumes that if inside mesh, radius is big enough to at least find 1 neighbor
        else {
            B.setZero(neighbors.size(), basis_dim);
            W.setZero(neighbors.size());
            F.setZero(neighbors.size());
            int k = 0;
            for (int i : neighbors) {
                Eigen::RowVector3d p_i = constrained_points.row(i);
                B.row(k) = get_basis_vector(p_i, basis_dim);
                W(k) = wendland((x - p_i).norm());
                F(k) = constrained_values(i);
                k++;
            }
            Eigen::MatrixXd M = W.asDiagonal() * B;
            Eigen::VectorXd V = W.asDiagonal() * F;
            Eigen::VectorXd c = M.fullPivHouseholderQr().solve(V); //TODO: fullPiv or colPiv?
            grid_values(i) = get_basis_vector(x, basis_dim).transpose().dot(c);
        }
    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "MLS computation time: " << duration.count() << endl;

    /*// Sphere center
    auto bb_min = grid_points.colwise().minCoeff().eval();
    auto bb_max = grid_points.colwise().maxCoeff().eval();
    Eigen::RowVector3d center = 0.5 * (bb_min + bb_max);

    double radius = 0.5 * (bb_max - bb_min).minCoeff();

    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);

    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x)
    {
        for (unsigned int y = 0; y < resolution; ++y)
        {
            for (unsigned int z = 0; z < resolution; ++z)
            {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);

                // Value at (x,y,z) = implicit function for the sphere
                grid_values[index] = (grid_points.row(index) - center).norm() - radius;
            }
        }
    }*/
}

void evaluateImplicitFunc_PolygonSoup()
{
    // Replace with your code here, for "key == '5'"
    //evaluateImplicitFunc();
    grid_values.resize(resolutionX * resolutionY * resolutionZ);
    Eigen::MatrixXd B;
    Eigen::VectorXd W;
    Eigen::VectorXd F;
    int basis_dim = get_basis_dim();
    for (int i = 0; i < grid_points.rows(); i++) {
        //Notations follow those from the course
        const auto x = grid_points.row(i);
        vector<int> neighbors = get_neighbors(x, true); //only original points
        if (neighbors.empty()) grid_values(i) = INT_MAX;
        else {
            B.setZero(neighbors.size(), basis_dim);
            W.setZero(neighbors.size());
            F.setZero(neighbors.size());
            int k = 0;
            for (int i : neighbors) {
                Eigen::RowVector3d p_i = constrained_points.row(i);
                B.row(k) = get_basis_vector(p_i, basis_dim);
                W(k) = wendland((x - p_i).norm());
                F(k) = S_k(x, i);
                k++;
            }
            Eigen::MatrixXd M = W.asDiagonal() * B;
            Eigen::VectorXd V = W.asDiagonal() * F;
            Eigen::VectorXd c = M.fullPivHouseholderQr().solve(V); //TODO: fullPiv or colPiv?
            grid_values(i) = get_basis_vector(x, basis_dim).transpose().dot(c);
        }
    }
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines()
{
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x < resolutionX; ++x)
    {
        for (unsigned int y = 0; y < resolutionY; ++y)
        {
            for (unsigned int z = 0; z < resolutionZ; ++z)
            {
                int index = x + resolutionX * (y + resolutionY * z);
                if (x < resolutionX - 1)
                {
                    int index1 = (x + 1) + y * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolutionY - 1)
                {
                    int index1 = x + (y + 1) * resolutionX + z * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolutionZ - 1)
                {
                    int index1 = x + y * resolutionX + (z + 1) * resolutionX * resolutionY;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

// Estimation of the normals via PCA.
void pcaNormal()
{
    //NP = -N; // to be replaced with your code
    NP.resize(N.rows(), 3);
    
    for(int i = 0; i < P.rows(); i++) {
        vector<int> neighbors = get_neighbors(P.row(i), true, pca_radius); //only original points
        Eigen::MatrixXd X_not_centered;
        X_not_centered.resize(neighbors.size(), 3);
        for(int j = 0; j < neighbors.size(); j++) {
            X_not_centered.row(j) = P.row(neighbors[j]);
        }
        Eigen::RowVector3d mu = X_not_centered.rowwise().mean();
        Eigen::MatrixXd X = X_not_centered.rowwise() - mu;
        Eigen::MatrixXd C = X.adjoint() * X;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(C);
        Eigen::MatrixXd eigenvectors = solver.eigenvectors();
        Eigen::RowVector3d normal = eigenvectors.col(0); //eigenvalues are sorted in increasing order
        double dot_prod = N.row(i).dot(normal);
        NP.row(i) = (dot_prod >= 0)? normal : -normal; //sign correction using provided normals
    }
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers)
{
    if (key == '1')
    {
        // Show imported points
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0, 0, 0));
    }

    if (key == '2')
    {
        // Show all constraints
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        // Add your code for computing auxiliary constraint points here
        constrained_points.setZero(3 * P.rows(), 3);
        constrained_values.setZero(3 * P.rows(), 1);
        double bounding_box_diag = (P.colwise().maxCoeff() - P.colwise().minCoeff()).norm();
        EPSILON = 0.01 * bounding_box_diag;
        create_spatial_index();
        auto start = high_resolution_clock::now();
        for (int i = 0; i < P.rows(); i++) {
            constrained_points.row(i) = P.row(i);
            constrained_values(i) = 0;
            Eigen::RowVector3d new_pos;
            double epsilon = EPSILON;
            double normal_norm = N.row(i).norm();
            auto normalized_normal = N.row(i) / normal_norm;
            do {
                new_pos = P.row(i) + epsilon * normalized_normal;
                epsilon /= 2;
            } while (!original_is_nearest(new_pos, i));
            constrained_points.row(i + P.rows()) = new_pos;
            constrained_values(i + P.rows()) = epsilon * 2;
            epsilon = EPSILON;
            do {
                new_pos = P.row(i) - epsilon * normalized_normal;
                epsilon /= 2;
            } while (!original_is_nearest(new_pos, i));
            constrained_points.row(i + 2 * P.rows()) = new_pos;
            constrained_values(i + 2 * P.rows()) = - epsilon * 2;
        }
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        cout << "constraints computation time: " << duration.count() << endl;
        // Add code for displaying all points, as above
        viewer.data().point_size = 7;
        viewer.data().add_points(constrained_points.block(0,            0,  P.rows(), 3), BLUE);
        viewer.data().add_points(constrained_points.block(P.rows(),     0,  P.rows(), 3), RED);
        viewer.data().add_points(constrained_points.block(2 * P.rows(), 0,  P.rows(), 3), GREEN);
    }

    if (key == '3')
    {
        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);
        // Add code for creating a grid
        // Add your code for evaluating the implicit function at the grid points
        // Add code for displaying points and lines
        // You can use the following example:

        /*** begin: sphere example, replace (at least partially) with your code ***/
        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc();

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i)
        {
            double value = grid_values(i);
            if (value < 0)
            {
                grid_colors(i, 1) = 1;
            }
            else
            {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }

        // Draw lines and points
        viewer.data().point_size = 6;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                                grid_lines.block(0, 3, grid_lines.rows(), 3),
                                Eigen::RowVector3d(0.8, 0.8, 0.8));
        /*** end: sphere example ***/
    }

    if (key == '4')
    {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0))
        {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolutionX, resolutionY, resolutionZ, V, F);
        if (V.rows() == 0)
        {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }
        Eigen::MatrixXd V_ = V;
        Eigen::MatrixXi F_ = F;
        if (connected_components) {
            run_connected_components(V_, F_);
        }

        igl::per_face_normals(V_, F_, FN);
        viewer.data().set_mesh(V_, F_);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);

        string output_file = "../res/" + file_name + ".off"; 
        igl::writeOFF(output_file, V_, F_);
        cout << "Wrote content to " << output_file << endl;
    }

    if (key == '5')
    {
        // Use the structure for key=='3' but replace the function evaluateImplicitFunc();
        // with a function performing the approximation of the implicit surface from polygon soup
        // Ref: Chen Shen, James F. O’Brien, and Jonathan Richard Shewchuk. Interpolating and approximating implicit surfaces from polygon soup.

        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core().align_camera_center(P);

        // Make grid
        createGrid();

        // Evaluate implicit function --> Function to be modified here
        evaluateImplicitFunc_PolygonSoup();

        // get grid lines
        getLines();

        // Display the reconstruction
        callback_key_down(viewer, '4', modifiers);
    }

    if (key == '6' || key == '7' || key == '8')
    {
        // Implement PCA Normal Estimation --> Function to be modified here
        pcaNormal();

        // To use the normals estimated via PCA instead of the input normals and then restaurate the input normals
        Eigen::MatrixXd N_tmp = N;
        N = NP;

        switch (key)
        {
        case '6':
            callback_key_down(viewer, '2', modifiers);
            break;
        case '7':
            callback_key_down(viewer, '3', modifiers);
            break;
        case '8':
            callback_key_down(viewer, '3', modifiers);
            callback_key_down(viewer, '4', modifiers);
            break;
        default:
            break;
        }

        // Restore input normals
        N = N_tmp;
    }

    return true;
}

bool callback_load_mesh(Viewer &viewer, string filename)
{
    igl::readOFF(filename, P, F, N);
    callback_key_down(viewer, '1', 0);
    return true;
}

int main(int argc, char *argv[])
{
    file_name = "sphere";
    if (argc != 2) {
        cout << "Usage ex2_bin <mesh_name (without path and without .off)>" << endl;
    } else {
        // Read points and normals
        file_name = argv[1];
    }
    igl::readOFF("../data/" + file_name + ".off", P, F, N);

    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        menu.draw_viewer_menu();

        // Add new group
        if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
        {
            // Expose variable directly ...
            //ImGui::InputInt("Resolution", &resolution, 0, 0);
            ImGui::InputInt("ResolutionX", &resolutionX, 0, 0);
            ImGui::InputInt("ResolutionY", &resolutionY, 0, 0);
            ImGui::InputInt("ResolutionZ", &resolutionZ, 0, 0);
            if (ImGui::Button("Reset Grid", ImVec2(-1, 0)))
            {
                std::cout << "ResetGrid\n";
                // Recreate the grid
                createGrid();
                // Switch view to show the grid
                callback_key_down(viewer, '3', 0);
            }

            // TODO: Add more parameters to tweak here...
            ImGui::InputDouble("EPSILON", &EPSILON);
            ImGui::SliderInt("polyDegree", &polyDegree, 0, 2, "%d", 0);
            ImGui::InputDouble("wendland radius", &wendlandRadius);
            ImGui::InputDouble("inflation ratio", &inflation_ratio);
            ImGui::Checkbox("PCA bounding box", &pca_bb);
            ImGui::InputDouble("PCA radius", &pca_radius);
            ImGui::Checkbox("Connected components", &connected_components);
            //TODO: Setup anisotropic grid resolution in the 3 axis
        }
    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}


/**
 * @brief 
 * Todo: In the report, insert comparison with/without connected components
 * 
 */