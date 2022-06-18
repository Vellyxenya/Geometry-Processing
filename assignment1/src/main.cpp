#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any libigl headers here ***/
#include <igl/vertex_triangle_adjacency.h>
#include <igl/per_corner_normals.h>
#include <igl/facet_components.h>
#include <igl/edge_topology.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Per-face normal array, #F x3
Eigen::MatrixXd FN;
// Per-vertex normal array, #V x3
Eigen::MatrixXd VN;
// Per-corner normal array, (3#F) x3
Eigen::MatrixXd CN;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd component_colors_per_face;

void subdivide_sqrt3(const Eigen::MatrixXd &V,
					 const Eigen::MatrixXi &F,
					 Eigen::MatrixXd &Vout,
					 Eigen::MatrixXi &Fout){
	// (1)
    Eigen::MatrixXd BC; // #F x 3
    int nb_original_rows = V.rows();
    std::vector<std::vector<int>> old_adj_vertices;
    igl::adjacency_list(F, old_adj_vertices);
    igl::barycenter(V, F, BC);
    Eigen::MatrixXd V_plus(V.rows() + BC.rows(), 3);
    V_plus << V, BC;
    Eigen::MatrixXi F3;
    F3.resize(F.rows() * 3, 3);
    for (int i = 0; i < F.rows(); i++) { //for each old face
        int v_center = nb_original_rows + i;
        int v0 = F(i, 0), v1 = F(i, 1), v2 = F(i, 2);
        F3.row(3*i  ) = Eigen::RowVector3i(v0, v1, v_center);
        F3.row(3*i+1) = Eigen::RowVector3i(v1, v2, v_center);
        F3.row(3*i+2) = Eigen::RowVector3i(v2, v0, v_center);
    }
    // (2)
    for (int i = 0; i < old_adj_vertices.size(); i++) {
        int n = old_adj_vertices[i].size();
        double a_n = (4 - 2 * std::cos(2 * M_PI / n)) / 9;
        Eigen::RowVector3d p = (1 - a_n) * V.row(i);
        Eigen::RowVector3d neighbor_vertices_sum(0, 0, 0);
        for (int j = 0; j < n; j++) {
            neighbor_vertices_sum += V.row(old_adj_vertices[i][j]);
        }
        p += neighbor_vertices_sum * a_n / n;
        V_plus.row(i) = p;
    }
    // (3) 
    Vout = V_plus;
    Fout = F3;
    Eigen::MatrixXi EV, FE, EF;
    igl::edge_topology(V_plus, F3, EV, FE, EF);
    for (int i = 0; i < EF.rows(); i++) {
        int face1 = EF(i, 0), face2 = EF(i, 1);
        if (face1 == -1 || face2 == -1) //border edge
            continue;
        int v1 = EV(i, 0), v2 = EV(i, 1);
        if (v1 >= nb_original_rows || v2 >= nb_original_rows) //new edge
            continue;
        //only old non-border edges remain
        Fout.row(face1) = Eigen::RowVector3i(F3(face1, 2), F3(face1, 0), F3(face2, 2));
        Fout.row(face2) = Eigen::RowVector3i(F3(face2, 2), F3(face2, 0), F3(face1, 2));
    }
}

//Custom variables
double cn_threshold = 50; //threshold for corner normals
bool print_vertices = true; //whether to print vertices when pressing '1' or '2'
Eigen::MatrixXd palette; //color palette for connected components

/**
 * @brief Utility function to print a 2d vector
 */
void print_v2(const std::vector<std::vector<int>>& v, const std::string type) {
    if (!print_vertices) return;
    for (int i = 0; i < v.size(); i++) {
        std::cout << type << i << ": ";
        for (int j = 0; j < v[i].size(); j++) {
            std::cout << v[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

/**
 * @brief Generate a random number between 0 and 1
 * @return double between 0 and 1
 */
double random_0_1() {
    return static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to face relations here;
        // store in VF,VFi.
        igl::vertex_triangle_adjacency(V, F, VF, VFi);
        print_v2(VF, "faces for v");
        std::cout << "----" << std::endl;
    }

    if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing vertex to vertex relations here:
        // store in VV.
        igl::adjacency_list(F, VV);
        print_v2(VV, "vertices for v");
        std::cout << "####" << std::endl;
    }

    if (key == '3') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        FN.setZero(F.rows(),3);
        // Add your code for computing per-face normals here: store in FN.
        igl::per_face_normals(V, F, FN);
        // Set the viewer normals.
        viewer.data().set_normals(FN);
    }

    if (key == '4') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-vertex normals here: store in VN.
        igl::per_vertex_normals(V, F, VN);
        // Set the viewer normals.
        viewer.data().set_normals(VN);
    }

    if (key == '5') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        // Add your code for computing per-corner normals here: store in CN.
        igl::per_corner_normals(V, F, cn_threshold, CN);
        //Set the viewer normals
        viewer.data().set_normals(CN);
    }

    if (key == '6') {
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        component_colors_per_face.setZero(F.rows(),3);
        // Add your code for computing per-face connected components here:
        // store the component labels in cid.
        igl::facet_components(F, cid);
        // Compute colors for the faces based on components, storing them in
        // component_colors_per_face.
        int nb_components = cid.maxCoeff() + 1;
        std::vector<int> nb_faces_per_component(nb_components, 0);
        palette.resize(nb_components, 3);
        for (int i = 0; i < nb_components; i++)
            palette.row(i) = Eigen::RowVector3d(random_0_1(), random_0_1(), random_0_1());
        for (int i = 0; i < cid.size(); i++) {
            component_colors_per_face.row(i) = palette.row(cid[i]);
            nb_faces_per_component[cid(i)]++;
        }
        for (int i = 0; i < nb_components; i++)
            std::cout << "Component " << (i+1) << " has " << nb_faces_per_component[i] << " faces." << std::endl;
        std::cout << "----" << std::endl;
        // Set the viewer colors
        viewer.data().set_colors(component_colors_per_face);
    }

    if (key == '7') {
		Eigen::MatrixXd Vout;
		Eigen::MatrixXi Fout;
        // Fill the subdivide_sqrt3() function with your code for sqrt(3) subdivision.
		subdivide_sqrt3(V, F, Vout, Fout);
        // Set up the viewer to display the new mesh
        V = Vout; F = Fout;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
    }

    return true;
}

bool load_mesh(Viewer& viewer,string filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
  igl::readOFF(filename,V,F);
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
        filename = std::string(argv[1]); // Mesh provided as command line argument
    }
    else {
        filename = std::string("../data/bunny.off"); // Default mesh
    }
	
    load_mesh(viewer,filename,V,F);

    callback_key_down(viewer, '1', 0);

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    // Standard code for drawing an additional window
    menu.callback_draw_custom_window = [&]() {
        ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiCond_FirstUseEver);
        ImGui::Begin(
            "Custom variables", nullptr,
            ImGuiWindowFlags_NoSavedSettings
        );
        ImGui::PushItemWidth(-80);
        ImGui::DragScalar("CN threshold", ImGuiDataType_Double, &cn_threshold, 0.1, 0, 0, "%.2f", 1.0F);
        ImGui::PopItemWidth();
        ImGui::Checkbox("print vertices", &print_vertices);
        ImGui::End();
    };
    
    viewer.launch();
}
