#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>

#include "Lasso.h"
#include "Colors.h"

//Custom includes
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>
#include <igl/cotmatrix.h>
#include <igl/slice.h>
#include <igl/local_basis.h>
#include <igl/grad.h>

//activate this for alternate UI (easier to debug)
//#define UPDATE_ONLY_ON_UP

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

//vertex array, #V x3
Eigen::MatrixXd V(0,3), V_cp(0, 3);
//face array, #F x3
Eigen::MatrixXi F(0,3);

//mouse interaction
enum MouseMode { SELECT, TRANSLATE, ROTATE, NONE };
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;

//for selecting vertices
std::unique_ptr<Lasso> lasso;
//list of currently selected vertices
Eigen::VectorXi selected_v(0,1);

//for saving constrained vertices
//vertex-to-handle index, #V x1 (-1 if vertex is free)
Eigen::VectorXi handle_id(0,1);
//list of all vertices belonging to handles, #HV x1
Eigen::VectorXi handle_vertices(0,1);
//centroids of handle regions, #H x1
Eigen::MatrixXd handle_centroids(0,3);
//updated positions of handle vertices, #HV x3
Eigen::MatrixXd handle_vertex_positions(0,3);
//index of handle being moved
int moving_handle = -1;
//rotation and translation for the handle being moved
Eigen::Vector3f translation(0,0,0);
Eigen::Vector4f rotation(0,0,0,1.);
typedef Eigen::Triplet<double> T;
//per vertex color array, #V x3
Eigen::MatrixXd vertex_colors;

//When false, use standard displacement vectors for details, when true use Deformation Transfer from part 2
bool use_deformation_transfer = false;

//function declarations (see below for implementation)
bool solve(Viewer& viewer);
void get_new_handle_locations();
Eigen::Vector3f computeTranslation (Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
Eigen::Vector4f computeRotation(Viewer& viewer, int mouse_x, int from_x, int mouse_y, int from_y, Eigen::RowVector3d pt3D);
void compute_handle_centroids();
Eigen::MatrixXd readMatrix(const char *filename);

bool callback_mouse_down(Viewer& viewer, int button, int modifier);
bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y);
bool callback_mouse_up(Viewer& viewer, int button, int modifier);
bool callback_pre_draw(Viewer& viewer);
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);
void onNewHandleID();
void applySelection();

//Custom definitions
//list of all vertices belonging to no handles (free vertices), #(V - HV) x1
Eigen::VectorXi V_f;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> solver_def_trans;

//Sliced matrices
Eigen::SparseMatrix<double> A_fc;
Eigen::SparseMatrix<double> B_f;
Eigen::SparseMatrix<double> M_fc;

Eigen::MatrixXd B, B_prime, S_prime; //The matrices designated by 'B', 'B'' and 'S'' in the assignment, respectively
Eigen::MatrixXd d_prime;
vector<int> longest_edges;
enum MeshMode { Mesh_B = 0, Mesh_B_prime, Mesh_S_prime, Mesh_S };
static MeshMode mesh_mode = Mesh_S_prime;
Eigen::SparseMatrix<double> GT_D;

void prefactorize() {
    //cout << "prefactorize" << endl;
    Eigen::SparseMatrix<double> L_w, M, A, M_inv, A_ff;
    igl::cotmatrix(V, F, L_w);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    igl::invert_diag(M, M_inv);
    A = L_w * M_inv * L_w;
    igl::slice(A, V_f, V_f, A_ff);
    igl::slice(A, V_f, handle_vertices, A_fc);

    solver.compute(A_ff); // the interior part of the (almost) bi-laplacian

    B = V;
    handle_vertex_positions = igl::slice(V, handle_vertices, 1);
    Eigen::MatrixXd b = -A_fc * handle_vertex_positions;
    Eigen::MatrixXd V_f_pos = solver.solve(b);
    igl::slice_into(V_f_pos, V_f, 1, B);

    Eigen::MatrixXd displacements = V - B;
    Eigen::MatrixXd X, Y, Z;
    igl::per_vertex_normals(B, F, Z);
    X.resize(V.rows(), 3);
    Y.resize(V.rows(), 3);
    vector<vector<int>> Adj;
    igl::adjacency_list(F, Adj);
    d_prime.resize(V.rows(), 3);
    longest_edges = vector<int>(V.rows());
    for (int i = 0; i < V.rows(); i++) {
        double max_squared_dist = 0;
        Eigen::RowVector3d x_i;
        int arg_neighbor = 0;
        for (int neighbor : Adj[i]) {
            double vertical_dist = (B.row(neighbor) - B.row(i)).dot(Z.row(i));
            Eigen::RowVector3d projected = B.row(neighbor) - vertical_dist * Z.row(i);
            Eigen::RowVector3d planar_vector = projected - B.row(i);
            double planar_squared_dist = planar_vector.squaredNorm();
            if (planar_squared_dist > max_squared_dist) {
                max_squared_dist = planar_squared_dist;
                x_i = planar_vector;
                arg_neighbor = neighbor;
            }
        }
        longest_edges[i] = arg_neighbor;
        x_i.normalize();
        Eigen::RowVector3d z_i = Z.row(i);
        Eigen::RowVector3d y_i = z_i.cross(x_i);
        double d_i_x = displacements.row(i).dot(x_i);
        double d_i_y = displacements.row(i).dot(y_i);
        double d_i_z = displacements.row(i).dot(z_i);
        d_prime.row(i) = Eigen::RowVector3d(d_i_x, d_i_y, d_i_z);
    }

    //### Deformation transfer
    //1) Compute the surface gradient
    Eigen::SparseMatrix<double> G;
    igl::grad(V, F, G);

    //2) Compute the weight matrix D
    Eigen::VectorXd double_area;
    igl::doublearea(V, F, double_area);
    Eigen::VectorXd replicated;
    replicated.resize(3 * F.rows());
    for(int i = 0; i < F.rows(); i++) {
      /*replicated(3*i) = double_area(i);
      replicated(3*i+1) = double_area(i);
      replicated(3*i+2) = double_area(i);*/
      replicated(i) = double_area(i);
      replicated(i+F.rows()) = double_area(i);
      replicated(i+2*F.rows()) = double_area(i);
    }
    //Note: D appears on both sides of the system so scaling does not matter
    Eigen::SparseMatrix<double> D = (Eigen::SparseMatrix<double>)(replicated.asDiagonal());

    //3) The Laplace and Divergence operators
    GT_D = G.transpose() * D;
    Eigen::SparseMatrix<double> GT_D_G = GT_D * G;

    Eigen::SparseMatrix<double> M_ff;
    igl::slice(GT_D_G, V_f, V_f, M_ff);
    igl::slice(GT_D_G, V_f, handle_vertices, M_fc);

    solver_def_trans.compute(M_ff);

    igl::slice(GT_D, V_f, 1, B_f);
}

bool solve(Viewer& viewer)
{
  Eigen::MatrixXd b = -A_fc * handle_vertex_positions;
  Eigen::MatrixXd V_f_pos = solver.solve(b);
  B_prime.resize(V.rows(), 3);
  igl::slice_into(V_f_pos, V_f, 1, B_prime);
  igl::slice_into(handle_vertex_positions, handle_vertices, 1, B_prime);

  if(!use_deformation_transfer) { //##### Multiresolution mesh editing
    Eigen::MatrixXd Z;
    igl::per_vertex_normals(B_prime, F, Z);
    S_prime.resize(V.rows(), 3);
    for (int i = 0; i < V.rows(); i++) {
        int j = longest_edges[i];
        Eigen::RowVector3d z = Z.row(i);
        double vertical_dist = (B_prime.row(j) - B_prime.row(i)).dot(z);
        Eigen::RowVector3d projected = B_prime.row(j) - vertical_dist * z;
        Eigen::RowVector3d x = (projected - B_prime.row(i)).normalized();
        Eigen::RowVector3d y = z.cross(x);
        Eigen::RowVector3d new_displacement = d_prime(i, 0) * x + d_prime(i, 1) * y + d_prime(i, 2) * z;
        S_prime.row(i) = B_prime.row(i) + new_displacement;
    }

    switch (mesh_mode) {
    case Mesh_B:
        V = B;
        break;
    case Mesh_B_prime:
        V = B_prime;
        break;
    case Mesh_S_prime:
        V = S_prime;
        break;
    case Mesh_S:
        V = V_cp;
    default:
        cout << "Not implemented" << endl;
        break;
    }
  } else { //##### Deformation transfer
    
    //1) Compute normals for 'B' and 'B''
    Eigen::MatrixXd B_prime_N;
    igl::per_face_normals(B_prime, F, B_prime_N);
    Eigen::MatrixXd B_N;
    igl::per_face_normals(B, F, B_N);

    //2) Compute S_j
    Eigen::MatrixX3d S_s(3 * F.rows(), 3);
    for(int i = 0; i < F.rows(); i++) {
      Eigen::Matrix3Xd s1;
      s1.resize(3, 3);
      s1.col(0) = (B_prime.row(F(i, 0)) - B_prime.row(F(i, 2))).transpose();
      s1.col(1) = (B_prime.row(F(i, 1)) - B_prime.row(F(i, 2))).transpose();
      s1.col(2) = (B_prime_N.row(i)).transpose();

      Eigen::Matrix3Xd s2;
      s2.resize(3, 3);
      s2.col(0) = (B.row(F(i, 0)) - B.row(F(i, 2))).transpose();
      s2.col(1) = (B.row(F(i, 1)) - B.row(F(i, 2))).transpose();
      s2.col(2) = (B_N.row(i)).transpose();

      Eigen::MatrixXd S_j = s1 * s2.inverse();
      Eigen::MatrixXd S_j_T = S_j.transpose();
      S_s.row(i) = S_j_T.row(0);
      S_s.row(i + F.rows()) = S_j_T.row(1);
      S_s.row(i + 2 * F.rows()) = S_j_T.row(2);
    }

    //4) Solve the linear system
    auto H = B_f * S_s;
    Eigen::MatrixX3d ps = solver_def_trans.solve(H - M_fc * handle_vertex_positions);

    //5) Update the vertices
    igl::slice_into(ps, V_f, 1, V);
    igl::slice_into(handle_vertex_positions, handle_vertices, 1, V);
  }

  return true;
};

void get_new_handle_locations()
{
  int count = 0;
  for (long vi = 0; vi < V.rows(); ++vi) {
    if (handle_id[vi] >= 0)
    {
      //cout << "if1" << endl;
      Eigen::RowVector3f goalPosition = V.row(vi).cast<float>();
      if (handle_id[vi] == moving_handle) {
        if (mouse_mode == TRANSLATE)
          goalPosition += translation;
        else if (mouse_mode == ROTATE) {
            //cout << "if2" << endl;
            /*if(moving_handle >= handle_centroids.rows())
              cout << "smash" << endl;*/
            Eigen::RowVector3f goalPositionCopy = goalPosition;
            goalPosition -= handle_centroids.row(moving_handle).cast<float>();
            //cout << goalPosition << rotation << goalPositionCopy << endl;
            //igl::rotate_by_quat(goalPosition.data(), rotation.data(), goalPositionCopy.data());
            Eigen::RowVector4f rotation_ = Eigen::RowVector4f(rotation);
            igl::rotate_by_quat(goalPosition.data(), rotation_.data(), goalPositionCopy.data());
            goalPosition = goalPositionCopy;
            goalPosition += handle_centroids.row(moving_handle).cast<float>();
            //cout << "if3" << endl;
        }
      }
      //cout << handle_vertex_positions.rows() << " " << count << endl;
      handle_vertex_positions.row(count++) = goalPosition.cast<double>();
      //cout << "if5" << endl;
    }
    //cout << "if6" << endl;
  }
  //cout << "after for" << endl;
}

bool load_mesh(string filename)
{
  igl::read_triangle_mesh(filename,V,F);
  viewer.data().clear();
  viewer.data().set_mesh(V, F);

  viewer.core().align_camera_center(V);
  V_cp = V;
  handle_id.setConstant(V.rows(), 1, -1);
  // Initialize selector
  lasso = std::unique_ptr<Lasso>(new Lasso(V, F, viewer));

  selected_v.resize(0,1);

  return true;
}

int main(int argc, char *argv[])
{
  if(argc != 2) {
    cout << "Usage assignment5 mesh.off>" << endl;
    load_mesh("../data/woody-lo.off");
  }
  else
  {
    load_mesh(argv[1]);
  }

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  menu.callback_draw_viewer_menu = [&]()
  {
    // Draw parent menu content
    menu.draw_viewer_menu();

    // Add new group
    if (ImGui::CollapsingHeader("Deformation Controls", ImGuiTreeNodeFlags_DefaultOpen))
    {
          int mouse_mode_type = static_cast<int>(mouse_mode);

          if (ImGui::Combo("Mouse Mode", &mouse_mode_type, "SELECT\0TRANSLATE\0ROTATE\0NONE\0"))
          {
            mouse_mode = static_cast<MouseMode>(mouse_mode_type);
          }

          if (ImGui::Button("Clear Selection", ImVec2(-1,0)))
          {
            selected_v.resize(0,1);
          }

          if (ImGui::Button("Apply Selection", ImVec2(-1,0)))
          {
            applySelection();
          }

          if (ImGui::Button("Clear Constraints", ImVec2(-1,0)))
          {
            handle_id.setConstant(V.rows(),1,-1);
          }
          if(ImGui::Checkbox("Deformation Transfer", &use_deformation_transfer)){}
    }

    if (ImGui::CollapsingHeader("Additional parameters", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Combo("Mesh Mode", (int*)(&mesh_mode), "Mesh_B\0Mesh_B_prime\0Mesh_S_prime\0Mesh_S\0\0");
    }
  };

  viewer.callback_key_down = callback_key_down;
  viewer.callback_mouse_down = callback_mouse_down;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_mouse_up = callback_mouse_up;
  viewer.callback_pre_draw = callback_pre_draw;

  viewer.data().point_size = 10;
  viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
  viewer.launch();
}


bool callback_mouse_down(Viewer& viewer, int button, int modifier)
{
  if (button == (int) Viewer::MouseButton::Right)
    return false;

  down_mouse_x = viewer.current_mouse_x;
  down_mouse_y = viewer.current_mouse_y;

  if (mouse_mode == SELECT)
  {
    if (lasso->strokeAdd(viewer.current_mouse_x, viewer.current_mouse_y) >=0)
      doit = true;
    else
      lasso->strokeReset();
  }
  else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
    //cout << "rotate" << endl;
    int vi = lasso->pickVertex(viewer.current_mouse_x, viewer.current_mouse_y);
    if(vi>=0 && handle_id[vi]>=0)  //if a region was found, mark it for translation/rotation
    {
      //cout << "rotate1" << endl;
      moving_handle = handle_id[vi];
      get_new_handle_locations();
      doit = true;
      //cout << "rotate2" << endl;
    }
  }
  return doit;
}

bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y)
{
  if (!doit)
    return false;
  if (mouse_mode == SELECT)
  {
    lasso->strokeAdd(mouse_x, mouse_y);
    return true;
  }
  if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
    if (mouse_mode == TRANSLATE) {
      translation = computeTranslation(viewer,
                                       mouse_x,
                                       down_mouse_x,
                                       mouse_y,
                                       down_mouse_y,
                                       handle_centroids.row(moving_handle));
    }
    else {
      rotation = computeRotation(viewer,
                                 mouse_x,
                                 down_mouse_x,
                                 mouse_y,
                                 down_mouse_y,
                                 handle_centroids.row(moving_handle));
    }
    get_new_handle_locations();
#ifndef UPDATE_ONLY_ON_UP
    solve(viewer);
    down_mouse_x = mouse_x;
    down_mouse_y = mouse_y;
#endif
    return true;

  }
  return false;
}

bool callback_mouse_up(Viewer& viewer, int button, int modifier)
{
  if (!doit)
    return false;
  doit = false;
  if (mouse_mode == SELECT)
  {
    selected_v.resize(0,1);
    lasso->strokeFinish(selected_v);
    return true;
  }

  if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
  {
#ifdef UPDATE_ONLY_ON_UP
    if(moving_handle>=0)
      solve(viewer);
#endif
    translation.setZero();
    rotation.setZero(); rotation[3] = 1.;
    moving_handle = -1;

    compute_handle_centroids();

    return true;
  }

  return false;
};


bool callback_pre_draw(Viewer& viewer)
{
  // initialize vertex colors
  vertex_colors = Eigen::MatrixXd::Constant(V.rows(),3,.9);

  // first, color constraints
  int num = handle_id.maxCoeff();
  if (num == 0)
    num = 1;
  for (int i = 0; i<V.rows(); ++i)
    if (handle_id[i]!=-1)
    {
      int r = handle_id[i] % MAXNUMREGIONS;
      vertex_colors.row(i) << regionColors[r][0], regionColors[r][1], regionColors[r][2];
    }
  // then, color selection
  for (int i = 0; i<selected_v.size(); ++i)
    vertex_colors.row(selected_v[i]) << 131./255, 131./255, 131./255.;

  viewer.data().set_colors(vertex_colors);
  viewer.data().V_material_specular.fill(0);
  viewer.data().V_material_specular.col(3).fill(1);
  viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE | igl::opengl::MeshGL::DIRTY_SPECULAR;


  //clear points and lines
  viewer.data().set_points(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXd::Zero(0,3));
  viewer.data().set_edges(Eigen::MatrixXd::Zero(0,3), Eigen::MatrixXi::Zero(0,3), Eigen::MatrixXd::Zero(0,3));

  //draw the stroke of the selection
  for (unsigned int i = 0; i<lasso->strokePoints.size(); ++i)
  {
    viewer.data().add_points(lasso->strokePoints[i],Eigen::RowVector3d(0.4,0.4,0.4));
    if(i>1)
      viewer.data().add_edges(lasso->strokePoints[i-1], lasso->strokePoints[i], Eigen::RowVector3d(0.7,0.7,0.7));
  }

  // update the vertex position all the time
  viewer.data().V.resize(V.rows(),3);
  viewer.data().V << V;

  viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_POSITION;

#ifdef UPDATE_ONLY_ON_UP
  //draw only the moving parts with a white line
  if (moving_handle>=0)
  {
    Eigen::MatrixXd edges(3*F.rows(),6);
    int num_edges = 0;
    for (int fi = 0; fi<F.rows(); ++fi)
    {
      int firstPickedVertex = -1;
      for(int vi = 0; vi<3 ; ++vi)
        if (handle_id[F(fi,vi)] == moving_handle)
        {
          firstPickedVertex = vi;
          break;
        }
      if(firstPickedVertex==-1)
        continue;


      Eigen::Matrix3d points;
      for(int vi = 0; vi<3; ++vi)
      {
        int vertex_id = F(fi,vi);
        if (handle_id[vertex_id] == moving_handle)
        {
          int index = -1;
          // if face is already constrained, find index in the constraints
          (handle_vertices.array()-vertex_id).cwiseAbs().minCoeff(&index);
          points.row(vi) = handle_vertex_positions.row(index);
        }
        else
          points.row(vi) =  V.row(vertex_id);

      }
      edges.row(num_edges++) << points.row(0), points.row(1);
      edges.row(num_edges++) << points.row(1), points.row(2);
      edges.row(num_edges++) << points.row(2), points.row(0);
    }
    edges.conservativeResize(num_edges, Eigen::NoChange);
    viewer.data().add_edges(edges.leftCols(3), edges.rightCols(3), Eigen::RowVector3d(0.9,0.9,0.9));

  }
#endif
  return false;

}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers)
{
  bool handled = false;
  if (key == 'S')
  {
    mouse_mode = SELECT;
    handled = true;
  }

  if ((key == 'T') && (modifiers == IGL_MOD_ALT))
  {
    mouse_mode = TRANSLATE;
    handled = true;
  }

  if ((key == 'R') && (modifiers == IGL_MOD_ALT))
  {
    mouse_mode = ROTATE;
    handled = true;
  }
  if (key == 'A')
  {
    applySelection();
    callback_key_down(viewer, '1', 0);
    handled = true;
  }
  // if (key == '2')
  // viewer.core().align_camera_position(B);
  //  }

  //viewer.ngui->refresh();
  return handled;
}

void onNewHandleID()
{
  //store handle vertices too
  int numFree = (handle_id.array() == -1).cast<int>().sum();
  int num_handle_vertices = V.rows() - numFree;
  handle_vertices.setZero(num_handle_vertices);
  V_f.setZero(numFree);
  handle_vertex_positions.setZero(num_handle_vertices,3);

  int count = 0, count_f = 0;
  for (long vi = 0; vi < V.rows(); ++vi)
    if (handle_id[vi] >= 0)
      handle_vertices[count++] = vi;
    else
      V_f[count_f++] = vi;

  compute_handle_centroids();
  prefactorize();
}

void applySelection()
{
  int index = handle_id.maxCoeff()+1;
  for (int i =0; i<selected_v.rows(); ++i)
  {
    const int selected_vertex = selected_v[i];
    if (handle_id[selected_vertex] == -1)
      handle_id[selected_vertex] = index;
  }
  selected_v.resize(0,1);

  onNewHandleID();
}

void compute_handle_centroids()
{
  //compute centroids of handles
  int num_handles = handle_id.maxCoeff()+1;
  handle_centroids.setZero(num_handles,3);

  Eigen::VectorXi num; num.setZero(num_handles,1);
  for (long vi = 0; vi<V.rows(); ++vi)
  {
    int r = handle_id[vi];
    if ( r!= -1)
    {
      handle_centroids.row(r) += V.row(vi);
      num[r]++;
    }
  }

  for (long i = 0; i<num_handles; ++i)
    handle_centroids.row(i) = handle_centroids.row(i).array()/num[i];

}

//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector3f computeTranslation (Viewer& viewer,
                                    int mouse_x,
                                    int from_x,
                                    int mouse_y,
                                    int from_y,
                                    Eigen::RowVector3d pt3D)
{
  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.data().model;
  //project the given point (typically the handle centroid) to get a screen space depth
  Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                      modelview,
                                      viewer.core().proj,
                                      viewer.core().viewport);
  float depth = proj[2];

  double x, y;
  Eigen::Vector3f pos1, pos0;

  //unproject from- and to- points
  x = mouse_x;
  y = viewer.core().viewport(3) - mouse_y;
  pos1 = igl::unproject(Eigen::Vector3f(x,y,depth),
                        modelview,
                        viewer.core().proj,
                        viewer.core().viewport);


  x = from_x;
  y = viewer.core().viewport(3) - from_y;
  pos0 = igl::unproject(Eigen::Vector3f(x,y,depth),
                        modelview,
                        viewer.core().proj,
                        viewer.core().viewport);

  //translation is the vector connecting the two
  Eigen::Vector3f translation = pos1 - pos0;
  return translation;

}


//computes translation for the vertices of the moving handle based on the mouse motion
Eigen::Vector4f computeRotation(Viewer& viewer,
                                int mouse_x,
                                int from_x,
                                int mouse_y,
                                int from_y,
                                Eigen::RowVector3d pt3D)
{

  Eigen::Vector4f rotation;
  rotation.setZero();
  rotation[3] = 1.;

  Eigen::Matrix4f modelview = viewer.core().view;// * viewer.data().model;

  //initialize a trackball around the handle that is being rotated
  //the trackball has (approximately) width w and height h
  double w = viewer.core().viewport[2]/8;
  double h = viewer.core().viewport[3]/8;

  //the mouse motion has to be expressed with respect to its center of mass
  //(i.e. it should approximately fall inside the region of the trackball)

  //project the given point on the handle(centroid)
  Eigen::Vector3f proj = igl::project(pt3D.transpose().cast<float>().eval(),
                                      modelview,
                                      viewer.core().proj,
                                      viewer.core().viewport);
  proj[1] = viewer.core().viewport[3] - proj[1];

  //express the mouse points w.r.t the centroid
  from_x -= proj[0]; mouse_x -= proj[0];
  from_y -= proj[1]; mouse_y -= proj[1];

  //shift so that the range is from 0-w and 0-h respectively (similarly to a standard viewport)
  from_x += w/2; mouse_x += w/2;
  from_y += h/2; mouse_y += h/2;

  //get rotation from trackball
  Eigen::Vector4f drot = viewer.core().trackball_angle.coeffs();
  Eigen::Vector4f drot_conj;
  igl::quat_conjugate(drot.data(), drot_conj.data());
  igl::trackball(w, h, float(1.), rotation.data(), from_x, from_y, mouse_x, mouse_y, rotation.data());

  //account for the modelview rotation: prerotate by modelview (place model back to the original
  //unrotated frame), postrotate by inverse modelview
  Eigen::Vector4f out;
  igl::quat_mult(rotation.data(), drot.data(), out.data());
  igl::quat_mult(drot_conj.data(), out.data(), rotation.data());
  return rotation;
}
