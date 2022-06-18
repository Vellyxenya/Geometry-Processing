#include <igl/read_triangle_mesh.h>
#include <igl/readTGF.h>
#include <igl/readDMAT.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/slice_into.h>
#include <igl/rotate_by_quat.h>
#include <igl/column_to_quats.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/jet.h>
#include <igl/dqs.h>
#include <igl/lbs_matrix.h>
#include <igl/harmonic.h>
#include <igl/quat_to_axis_angle.h>
#include <igl/polar_dec.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/project_to_line_segment.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>
#include <igl/cotmatrix.h>
#include <igl/slice.h>
#include <igl/local_basis.h>
#include <igl/grad.h>

#include "Colors.h"

#include <math.h>
#include <imgui/imgui.h>
#include <Eigen/IterativeLinearSolvers>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;
using vecQ = std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>>;

//Colors
#define RED 1.0, 0, 0
#define PURPLE 1.0, 0.1, 1.0
#define GREENISH 0.3, 1.0, 0.6
#define WHITISH 0.9, 0.9, 0.9

Viewer viewer;
string folder; //folder to load the data from

Eigen::MatrixXd V(0, 3); //vertex array, #V x3
Eigen::MatrixXi F(0, 3); //face array, #F x3
Eigen::MatrixXd vertex_colors; //per vertex color array, #V x3
Eigen::MatrixXd Joints__; //Rest Joint positions
Eigen::MatrixXd Joints_; //Current Joint positions
Eigen::MatrixXi Links; //Links (Edges between joints)
Eigen::MatrixXd Matrot_; //Stacked Rotation matrices (per frame and then per joint)
Eigen::MatrixXi Parents; //Joint parent indices
Eigen::VectorXi H_; //Handle indices
Eigen::MatrixXd HarmonicWeights_; //vertex harmonic weights
Eigen::MatrixXd FaceWeights_; //face weights computed from vertex weights
Eigen::MatrixXd V_; //currently displayed vertices
Eigen::MatrixXd V_root; //vertices assigned to the root handle
Eigen::VectorXi root_vertices, non_root_vertices; //indices of root/non-root vertices
Eigen::VectorXi handle_vertices; //non-free vertices
Eigen::VectorXi free_vertices; //free vertices
Eigen::SparseMatrix<double> G_; //Gradient matrix
Eigen::SparseMatrix<double> G_c; //constrained elements of the gradient matrix
Eigen::RowVector3d Root; //Root position of base mesh
Eigen::SparseMatrix<double> M_fc; //Submatrix used in poisson-stitching
Eigen::SparseMatrix<double> B_f; //Submatrix used in poisson-stitching
Eigen::MatrixXd joint_color(1, 3);
Eigen::MatrixXd link_color(1, 3);

//vectors
vecQ vQuat_;
vecQ vIdQuat;
vector<vecQ> example_poses;
vector<Eigen::Vector3d> vTrans_;
vector<Eigen::MatrixXd> Displacements_;
vector<Eigen::MatrixXd> Deformations_;
vector<Eigen::MatrixXd> unposed_examples_per_vertex;
vector<Eigen::MatrixXd> unposed_examples_per_face;

//enums
enum SkinningMode { PER_VERTEX, DUAL_QUATERNION, PER_FACE };
SkinningMode skinning_mode = PER_VERTEX;
enum RBF_Function { GAUSSIAN, POLYHARMONIC, THIN_PLATE };
RBF_Function rbf_function = POLYHARMONIC;
enum HandlesColorMode { ALL_HANDLES, HARMONIC_SINGLE_HANDLE, SINGLE_HANDLE };
HandlesColorMode handles_color_mode = ALL_HANDLES;

//bools
bool use_quaternions = false; //Use the provided quaternions if available or not
bool animate = false; //Animate or pause the animation
bool use_reference_handles = true; //whether to use reference handles or not
bool example_data_available = true; //Data files are available
bool quaternions_available = false; //Quaternions file is available
bool context_aware = false; //Context-aware animation
bool visualize_pose = false; //whether to visualize an example pose instead of animating
bool equation_10 = false; //whether to use equation 9 or 10 for per-face context-aware
bool bidirectional = false; //whether to animate backwards when end of animation reached
bool forward_animation = true; //animate forwards of backwards
bool per_vertex_unposing = true; //per-vertex or per-face based unposing

int K, L; //Nb Joints and nb frames
int J = 0; //Nb of examples;
int l = 0; //frame counter
int selected_handle = 0; //currently highlighted handle
int root_link;
int current_pose = 0; //example pose displayed on screen
double epsilon = 5; //Gaussian RBF hyper param

//Solvers
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> laplaceSolver;
Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>, Eigen::RowMajor> faceLBSSolver;

void draw_viewer_menu_minimal(Viewer* viewer);
bool callback_mouse_down(Viewer& viewer, int button, int modifier);
bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y);
bool callback_mouse_up(Viewer& viewer, int button, int modifier);
bool callback_pre_draw(Viewer& viewer);
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);
void compute_harmonic_weights();
void per_face_lbs();
void blend_quaternions(const vecQ& quats, const Eigen::VectorXd& weights,
  Eigen::Quaterniond& quat);
void blend_translations(const vector<Eigen::Vector3d>& translations,
  const Eigen::VectorXd& weights, Eigen::Vector3d& translation);


//RBF kernels for context-based deformation
double RBF(double r) {
  switch (rbf_function) {
  case GAUSSIAN:
    return exp(-pow(epsilon * r, 2));
    break;
  case POLYHARMONIC:
    return pow(r, 3);
    break;
  case THIN_PLATE:
    return r == 0 ? 0 : r * r * log(r);
    break;
  default:
    cout << "Not implemented" << endl;
    return 0;
    break;
  }
}

//cf.: https://math.stackexchange.com/questions/90081/quaternion-distance
double quaternion_distance(const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2) {
  return 1 - pow(q1.dot(q2), 2); //rough distance between quaternions
}

//cumulated distance between 2 lists of quaternions (2 poses)
double quaternion_distance(const vecQ& q1, const vecQ& q2) {
  double distance = 0;
  for(int i = 0; i < q1.size(); i++) {
    distance += pow(quaternion_distance(q1[i], q2[i]), 2);
  }
  return sqrt(distance);
}

//Given a pose, compute the weights representing how close the pose is
//to reference poses
void compute_example_weights(const vecQ& pose, 
  const vector<vecQ>& example_poses, Eigen::RowVectorXd& a_js) {
  //Compute Kernel matrix
  a_js.resize(J);
  Eigen::MatrixXd Phi(J, J);
  for(int j = 0; j < J; j++) {
    for(int t = 0; t < J; t++) {
      Phi(j, t) = RBF(quaternion_distance(example_poses[j], example_poses[t]));
    }
  }
  Eigen::LLT<Eigen::MatrixXd> phiSolver;
  phiSolver.compute(Phi);
  double sum = 0;
  for(int j = 0; j < J; j++) { //For each example, we solve a linear system to compute a_j(P)
    //Set the constraints
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(J);
    rhs(j) = 1;
    //solve for the c_js
    Eigen::VectorXd c_j = phiSolver.solve(rhs);
    //compute the corresponding weights
    double a_j = 0;
    for(int t = 0; t < J; t++) {
      a_j += c_j(t) * RBF(quaternion_distance(pose, example_poses[t]));
    }
    a_js(j) = a_j;
    sum += a_j;
  }
  //Normalize
  for(int j = 0; j < J; j++) {
    a_js(j) /= sum;
  }
}

//Unpose the vertices of a given mesh given absolute rotations and translations
//handles both vertex-based and face-based unposing
void unpose(const Eigen::MatrixXd& example_mesh, const vecQ& quats,
  const std::vector<Eigen::Vector3d>& translations,
  Eigen::MatrixXd& unposed_per_vertex, Eigen::MatrixXd& unposed_per_face) {
  //Per-vertex unposing
  unposed_per_vertex.setZero(V.rows(), 3);
  for(int i = 0; i < V.rows(); i++) {
    for(int j = 0; j < K; j++) {
      unposed_per_vertex.row(i) += HarmonicWeights_(i, j) * (quats[j].inverse() * (example_mesh.row(i) - translations[j].transpose()));
    }
  }
  //Per-face unposing
  unposed_per_face.setZero(V.rows(), 3);
  Eigen::MatrixXd Q_Tilde(3 * F.rows(), 3);
  for(int i = 0; i < F.rows(); i++) {
    Eigen::RowVector3d v0 = Eigen::RowVector3d::Zero();
    Eigen::RowVector3d v1 = Eigen::RowVector3d::Zero();
    Eigen::RowVector3d v2 = Eigen::RowVector3d::Zero();
    for(int j = 0; j < K; j++) { //unpose each vertex of the face
      v0 += FaceWeights_(i, j) * (quats[j].inverse() * (example_mesh.row(F(i, 0)) - translations[j].transpose()));
      v1 += FaceWeights_(i, j) * (quats[j].inverse() * (example_mesh.row(F(i, 1)) - translations[j].transpose()));
      v2 += FaceWeights_(i, j) * (quats[j].inverse() * (example_mesh.row(F(i, 2)) - translations[j].transpose()));
    }

    Eigen::Matrix3Xd s1(3, 3);
    s1.col(0) = (v0 - v2).transpose();
    s1.col(1) = (v1 - v2).transpose();
    s1.col(2) = s1.col(0).cross(s1.col(1)).normalized();

    Eigen::Matrix3Xd s2(3, 3);
    s2.col(0) = (V.row(F(i, 0)) - V.row(F(i, 2))).transpose();
    s2.col(1) = (V.row(F(i, 1)) - V.row(F(i, 2))).transpose();
    s2.col(2) = s2.col(0).cross(s2.col(1)).normalized();

    //Compute 3x3 transformation matrices
    Eigen::MatrixXd q_tilde_j = s1 * s2.inverse();
    Eigen::MatrixXd q_tilde_j_T = q_tilde_j.transpose();
    Q_Tilde.row(i) = q_tilde_j_T.row(0);
    Q_Tilde.row(i + F.rows()) = q_tilde_j_T.row(1);
    Q_Tilde.row(i + 2 * F.rows()) = q_tilde_j_T.row(2);
  }
  //Solve the Poisson-stitching linear system
  Eigen::MatrixXd V_non_root = faceLBSSolver.solve(B_f * Q_Tilde - M_fc * V_root);
  igl::slice_into(V_root, root_vertices, 1, unposed_per_face);
  igl::slice_into(V_non_root, non_root_vertices, 1, unposed_per_face);
}

//Transform a list of rotation matrices to quaternions
vecQ quats_from_rots(Eigen::MatrixXd rots) {
  vecQ pose;
  for(int k = 0; k < K; k++) {
    Eigen::Matrix3d rot_mat = rots.block(3*k, 0, 3, 3);
    Eigen::Quaterniond quat = Eigen::Quaterniond(rot_mat);
    pose.push_back(quat);
  }
  return pose;
}

void load_example_data(vector<vecQ>& example_poses) {
  cout << "Loading example data" << endl;
  Eigen::MatrixXd junk;
  J = 0; //example index
  while(true) { //while examples files are available
    Eigen::MatrixXd Matrot_eg;
    //If reading fails, there are no more examples
    if(!igl::readDMAT("../data/"+folder+"/eg"+std::to_string(J)+".dmat", Matrot_eg)) return;
    Eigen::MatrixXd example_mesh;
    if(!igl::read_triangle_mesh("../data/"+folder+"/eg"+std::to_string(J)+".obj", example_mesh, junk)) return;
    J++;
    //Center the examples like we did for the original mesh
    example_mesh.rowwise() -= Root;
    vecQ pose = quats_from_rots(Matrot_eg);
    vecQ vQuats;
    std::vector<Eigen::Vector3d> vTrans;
    //Perform forward kinematics to compute the absolute rotations and translations
    igl::forward_kinematics(Joints__, Links, Parents, pose, vTrans_, vQuats, vTrans);
    example_poses.push_back(pose);
    Eigen::MatrixXd V_unposed_per_vertex, V_unposed_per_face;
    //Unpose using computed rotations and trnaslations
    unpose(example_mesh, vQuats, vTrans, V_unposed_per_vertex, V_unposed_per_face);
    unposed_examples_per_vertex.push_back(V_unposed_per_vertex);
    unposed_examples_per_face.push_back(V_unposed_per_face);

    //Per-vertex: We simply need the displacements
    Displacements_.push_back(V_unposed_per_vertex - V);

    //Per-face: we need to compute the deformations for Poisson stitching
    Eigen::MatrixXd Q_Tilde(3 * F.rows(), 3);
    for(int i = 0; i < F.rows(); i++) {
      Eigen::Matrix3Xd s1(3, 3);
      s1.col(0) = (V_unposed_per_face.row(F(i, 0)) - V_unposed_per_face.row(F(i, 2))).transpose();
      s1.col(1) = (V_unposed_per_face.row(F(i, 1)) - V_unposed_per_face.row(F(i, 2))).transpose();
      s1.col(2) = s1.col(0).cross(s1.col(1)).normalized();

      Eigen::Matrix3Xd s2(3, 3);
      s2.col(0) = (V.row(F(i, 0)) - V.row(F(i, 2))).transpose();
      s2.col(1) = (V.row(F(i, 1)) - V.row(F(i, 2))).transpose();
      s2.col(2) = s2.col(0).cross(s2.col(1)).normalized();

      Eigen::MatrixXd q_tilde_j = s1 * s2.inverse();
      //Careful with indexing here, we keep the matrices contiguous
      Q_Tilde.block(3*i, 0, 3, 3) = q_tilde_j;
    }
    Deformations_.push_back(Q_Tilde);
  }
}

//Compute distance between some vertex to some bone if it's inside a cylinder
//around the bone
double distance_vertex_to_link(const Vector3d& v, const Vector2i& link, double influence = 1) {
  Vector3d v1 = Joints__.row(link(0)), v2 = Joints__.row(link(1));
  //create a cylinder that is smaller or equal in length to the cylinder defined by the link
  //this is to reduce self-collisions at the joints
  v2 = v1 + (v2 - v1) * influence;
  v1 = v2 + (v1 - v2) * influence;
  double link_length = (v2 - v1).squaredNorm();
  double t = (v - v1).dot(v2 - v1) / link_length;
  if(t < 0 || t > 1) return 100000; //Falls ouside the cylinder so return large distance
  Eigen::Vector3d projection = v1 + t * (v2 - v1); // Projection falls onto the segment
  return (v - projection).norm();
}

//Either load reference handles or compute ones from scratch
void compute_handles(Eigen::VectorXi& H_) {
  if(use_reference_handles) {
    cout << "using reference handles" << endl;
    igl::readDMAT("../data/"+folder+"/handles.dmat", H_);
  } else {
    cout << "using computed handles" << endl;
    H_.resize(V.rows());
    vector<double> min_dists;
    //Compute minimum distance between links and nearest vertex
    for(int k = 0; k < K; k++) {
      double min_dist_k = 1000000;
      for(int i = 0; i < V.rows(); i++) {
        double dist = distance_vertex_to_link(V.row(i), Links.row(k));
        if(dist < min_dist_k)
          min_dist_k = dist;
      }
      min_dists.push_back(min_dist_k);
    }
    //For each vertex, assign it to the closest handle if it is close enough
    for(int i = 0; i < V.rows(); i++) {
      double min_dist = 1000000;
      int argmin_k = -1;
      for(int k = 0; k < K; k++) {
        double dist = distance_vertex_to_link(V.row(i), Links.row(k), 0.9);
        //Only consider vertex if it falls within some threshold
        //The factor 2.1 is empirical and seems to work best
        if(dist < min_dists[k] * 2.1 && dist < min_dist) {
          min_dist = dist;
          argmin_k = k;
        }
      }
      H_(i) = argmin_k;
    }
    //### Post-Processing ###
    //If 2 adjacent handles touch, set the vertex in question to be free
    //This reduces self-intersections between vertices from different handles
    vector<vector<int>> A;
    igl::adjacency_list(F, A);
    for(int i = 0; i < A.size(); i++) {
      int handle = H_(i);
      if(handle == -1) continue;
      for(int j = 0; j < A[i].size(); j++) {
        if(H_(A[i][j]) != -1 && H_(A[i][j]) != handle) {
          H_(i) = -1;
          H_(A[i][j]) = -1;
        }
      }
    }
    //We create even further distance between regions: if some free vertex
    //has at least 2 neighboring vertices with different handles, set those
    //said vertices to free.
    for(int i = 0; i < A.size(); i++) {
      int handle = H_(i);
      if(handle != -1) continue;
      std::set<int> handles;
      for(int j = 0; j < A[i].size(); j++) {
        if(H_(A[i][j]) != -1) {
          handles.insert(H_(A[i][j]));
        }
      }
      if(handles.size() > 1) {
        for(int j = 0; j < A[i].size(); j++) {
          H_(A[i][j]) = -1;
        }
      }
    }
  }

  //Compute free vertices, handles vertices, root vertices etc.
  int num_free = (H_.array() == -1).cast<int>().sum();
  int num_handle_vertices = H_.size() - num_free;
  int num_root_vertices = (H_.array() == root_link).cast<int>().sum();
  int num_non_root_vertices = H_.size() - num_root_vertices;
  handle_vertices.setZero(num_handle_vertices);
  free_vertices.setZero(num_free);
  root_vertices.setZero(num_root_vertices);
  non_root_vertices.setZero(num_non_root_vertices);
  int count = 0, count_f = 0, count_r = 0, count_nr = 0;
  for (long vi = 0; vi < H_.size(); ++vi) {
    if (H_(vi) >= 0) {
      handle_vertices[count++] = vi;
      if(H_(vi) == root_link)
        root_vertices[count_r++] = vi;
      else
        non_root_vertices[count_nr++] = vi;
    } else {
      free_vertices[count_f++] = vi;
      non_root_vertices[count_nr++] = vi;
    }
  }
}

//This processing depends on the assigned handles and should hence
//Be re-executed each time the handles are changed
void handle_dependent_processing() {
  compute_handles(H_);
  compute_harmonic_weights();
  per_face_lbs();
  load_example_data(example_poses);
  example_data_available = J > 0;
}

//Load files and perform some preprocessing
bool load_data() {
  igl::readTGF("../data/"+folder+"/skeleton.tgf", Joints__, Links);
  Root = Joints__.row(0).eval();
  //Center the joints as it makes things easier
  for(int i = 0; i < Joints__.rows(); i++) {
    Joints__.row(i) = Joints__.row(i).eval() - Root;
  }
  Joints_ = Joints__; //Current Joint positions
  igl::directed_edge_parents(Links, Parents);
  //Determine the root link
  for(int i = 0; i < Links.rows(); i++) {
    if(Parents(i) == -1) {
      root_link = i;
      break;
    }
  }
  K = Links.rows(); //Nb joints
  
  //Compute the translation vector to be used in forward kinematics
  //Just a zero vector actually
  vTrans_ = vector<Eigen::Vector3d>(Links.rows(), Eigen::Vector3d::Zero());

  //The actual mesh
  if(igl::read_triangle_mesh("../data/"+folder+"/reference.off", V, F)) {}
  else igl::read_triangle_mesh("../data/"+folder+"/reference.obj", V, F);
  V.rowwise() -= Root; //Offset by the root as it makes everything easier
  viewer.data().clear();
  viewer.data().set_mesh(V, F);
  viewer.core().align_camera_center(V);
  V_ = Eigen::MatrixXd(V.rows(), 3);

  igl::readDMAT("../data/"+folder+"/rots_all_frames.dmat", Matrot_);
  L = Matrot_.rows()/ (3 * K); //Nb frames

  Eigen::VectorXd Quat_; //If we have a quaternion data file, load it
  if(quaternions_available = igl::readDMAT("../data/"+folder+"/pose_quats.dmat", Quat_)) {
    igl::column_to_quats(Quat_, vQuat_);
    vIdQuat = std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>>(K, Eigen::Quaterniond::Identity());
  }

  //Computing vertex harmonic weights, face weights, unposing examples
  //etc. is handle dependent
  handle_dependent_processing();

  //Sanity check: if we set the pose to say example_poses[2], we should get
  //the weights [0, 0, 1, 0]
  //Careful, for the Gaussian Kernel, setting epsilon too small
  //messes things up for some reason!
  //Just increase epsilon until the sanity check passes
  /*if(example_data_available) {
    vector<double> example_weights;
    compute_example_weights(example_poses[2], example_poses, example_weights);
    for(int i = 0; i < example_weights.size(); i++) {
      cout << "weight: " << example_weights[i] << endl;
    }
  }*/

  return true;
}

//Adaped from igl::forward_kinematics
void forward_kinematics(const Eigen::MatrixXd & C, const Eigen::MatrixXi & BE, const Eigen::VectorXi & P,
  const std::vector<Eigen::Matrix3d> & dQ, const std::vector<Eigen::Vector3d> & dT,
  std::vector<Eigen::Matrix3d> & vQ, std::vector<Eigen::Vector3d> & vT) {
  const int m = BE.rows(); 
  assert(m == P.rows());
  assert(m == (int)dQ.size());
  assert(m == (int)dT.size());
  vector<bool> computed(m, false);
  vQ.resize(m);
  vT.resize(m);
  // Dynamic programming
  function<void (int)> fk_helper = [&] (int b) {
    if(!computed[b]) {
      if(P(b) < 0) {
        // base case for roots
        vQ[b] = dQ[b];
        const Vector3d r = C.row(BE(b,0)).transpose();
        vT[b] = r-dQ[b]*r + dT[b];
      } else {
        // Otherwise first compute parent's
        const int p = P(b);
        fk_helper(p);
        vQ[b] = vQ[p] * dQ[b];
        const Vector3d r = C.row(BE(b,0)).transpose();
        vT[b] = vT[p] - vQ[b]*r + vQ[p]*(r + dT[b]);
      }
      computed[b] = true;
    }
  };
  for(int b = 0; b<m; b++) {
    fk_helper(b);
  }
}

//Standard quaternion slerp
Eigen::Quaterniond slerp(const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2, double t) {
  return q1.slerp(t, q2);
}

//Slerp between 2 poses consisting in lists of quaternions
vecQ slerp(vecQ& q1, vecQ& q2, double t) {
  vecQ v;
  v.resize(q1.size());
  for(int i = 0; i < v.size(); i++) { //Slerp each pair of quaternions
    v[i] = slerp(q1[i], q2[i], t);
  }
  return v;
}

//Compute per-vertex harmonic weights
void compute_harmonic_weights() {
  //Compute the Laplacian and only use its free part
  Eigen::SparseMatrix<double> L, L_ff, L_fc;
  igl::cotmatrix(V, F, L);
  igl::slice(L, free_vertices, free_vertices, L_ff);
  igl::slice(L, free_vertices, handle_vertices, L_fc);
  laplaceSolver.compute(L_ff);

  //Initialize the #V x #K weights matrix (each vertex has K weights)
  Eigen::SparseMatrix<double> WSparse;
  Eigen::MatrixXd Wc;
  WSparse.resize(H_.size(), K);
  vector<Eigen::Triplet<double>> data;
  for(int i = 0; i < H_.size(); i++) { //for each vertex
    if(H_(i) < 0) continue; //free vertex so do nothing
    data.push_back(Eigen::Triplet<double>(i, H_(i), 1));
  }
  WSparse.setFromTriplets(data.begin(), data.end());
  Eigen::MatrixXd W = WSparse.toDense();
  igl::slice(W, handle_vertices, 1, Wc);

  //Solve for the weights of the free vertices
  Eigen::MatrixXd Wf = laplaceSolver.solve(-L_fc * Wc); //fxK matrix, where f is # free vertices

  //Remerge the free and non-free weights in a single matrix
  HarmonicWeights_.resize(V.rows(), K);
  igl::slice_into(Wf, free_vertices, 1, HarmonicWeights_);
  igl::slice_into(Wc, handle_vertices, 1, HarmonicWeights_);

  //Sanity check: indeed the sum over the rows gives 1 for each row
  //cout << HarmonicWeights_.rowwise().sum().transpose() << endl;
}

//Precompute matrices to solve for the per-face linear blend skinning 
void per_face_lbs() {
  //Compute face weights given vertex weights
  FaceWeights_.resize(F.rows(), K);
  for(int i = 0; i < F.rows(); i++) {
    auto row = F.row(i);
    int v0 = row(0), v1 = row(1), v2 = row(2);
    FaceWeights_.row(i) = (HarmonicWeights_.row(v0) + HarmonicWeights_.row(v1) + HarmonicWeights_.row(v2)) / 3;
  }

  //Compute the gradient matrix of the mesh
  igl::grad(V, F, G_);

  //Compute the weight matrix D
  Eigen::VectorXd double_area;
  igl::doublearea(V, F, double_area);
  Eigen::VectorXd replicated;
  replicated.resize(3 * F.rows());
  for(int i = 0; i < F.rows(); i++) {
    replicated(i) = double_area(i);
    replicated(i+F.rows()) = double_area(i);
    replicated(i+2*F.rows()) = double_area(i);
  }
  //Note: D appears on both sides of the system so scaling does not matter
  Eigen::SparseMatrix<double> D = (Eigen::SparseMatrix<double>)(replicated.asDiagonal());

  //Compute the Laplace and Divergence operators
  Eigen::SparseMatrix<double> GT_D = G_.transpose() * D;
  Eigen::SparseMatrix<double> GT_D_G = GT_D * G_;

  //Factorize the linear system
  Eigen::SparseMatrix<double> M_ff;
  igl::slice(GT_D_G, non_root_vertices, non_root_vertices, M_ff);
  igl::slice(GT_D_G, non_root_vertices, root_vertices, M_fc);
  igl::slice(V, root_vertices, 1, V_root);
  igl::slice(GT_D, non_root_vertices, 1, B_f);
  igl::slice(V, root_vertices, 1, V_root);
  faceLBSSolver.analyzePattern(M_ff);
  faceLBSSolver.compute(M_ff);
}

//Compute the natural logarithm of a quaternion
//Assumes a unit quaternion input
Eigen::Quaterniond compute_log_quat(const Eigen::Quaterniond& quat) {
  Eigen::Vector3d vec = quat.vec().normalized() * acos(quat.w());
  return Eigen::Quaterniond(0, vec(0), vec(1), vec(2));
}

//Compute the exponential of a quaternion
//Assumes a unit quaternion input
Eigen::Quaterniond compute_exp_quat(const Eigen::Quaterniond& quat) {
  double v_norm;
  if((v_norm = quat.vec().norm()) == 0) return Eigen::Quaterniond(exp(1), 0, 0, 0);
  double exp_a = exp(quat.w());
  Eigen::Vector3d vec = exp_a * sin(v_norm) * quat.vec().normalized();
  return Eigen::Quaterniond(exp_a * cos(v_norm), vec(0), vec(1), vec(2));
}

//blend quaternions in log-space using link weights
void blend_quaternions(const vecQ& quats, const Eigen::VectorXd& weights, Eigen::Quaterniond& quat) {
  quat = Eigen::Quaterniond(0, 0, 0, 0);
  for(int k = 0; k < weights.size(); k++) {
    quat.coeffs() += weights(k) * compute_log_quat(quats[k]).coeffs();
  }
  quat = compute_exp_quat(quat);
  quat.normalize();
}

//blend translations given link weights
void blend_translations(const vector<Eigen::Vector3d>& translations,
  const Eigen::VectorXd& weights, Eigen::Vector3d& translation) {
  translation = Eigen::Vector3d::Zero();
  for(int k = 0; k < K; k++) {
    translation += weights(k) * translations[k];
  }
}

void highlight_handles() {
  vertex_colors = Eigen::MatrixXd::Constant(V.rows(), 3, 0.9);
  switch (handles_color_mode) {
  case ALL_HANDLES: //Highlight all handles (no skinning weights)
    for (int i = 0; i < V.rows(); ++i) {
      if (H_(i) != -1) {
        int color_id = H_(i) % MAXNUMREGIONS;
        vertex_colors.row(i) << regionColors[color_id][0], regionColors[color_id][1], regionColors[color_id][2];
      }
    }
    break;
  case HARMONIC_SINGLE_HANDLE: //Only highlight a single handle (with skinning weights)
    for (int i = 0; i < V.rows(); ++i) {
      double w_k = HarmonicWeights_(i, selected_handle);
      int color_id = selected_handle % MAXNUMREGIONS;
      vertex_colors.row(i) << regionColors[color_id][0], regionColors[color_id][1], regionColors[color_id][2];
      vertex_colors.row(i) *= w_k;
      vertex_colors.row(i) += (1 - w_k) * Eigen::RowVector3d(WHITISH);
    }
    break;
  case SINGLE_HANDLE: //Only highlight a single handle (no skinning weights)
    for (int i = 0; i < V.rows(); ++i) {
      double w_k = HarmonicWeights_(i, selected_handle);
      int color_id = selected_handle % MAXNUMREGIONS;
      if(w_k == 1)
        vertex_colors.row(i) << regionColors[color_id][0], regionColors[color_id][1], regionColors[color_id][2];
    }
    break;
  default:
    cout << "Handle color mode not implemented" << endl;
    break;
  }
  
  viewer.data().set_colors(vertex_colors);
  viewer.data().V_material_specular.fill(0);
  viewer.data().V_material_specular.col(3).fill(1);
  viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE | igl::opengl::MeshGL::DIRTY_SPECULAR;
}

//Update the underlying rig (joint positions) using forward kinematics
void update_rig(vecQ& lQuats, vecQ& vQuats, std::vector<Eigen::Vector3d>& vTrans,
  std::vector<Eigen::Matrix3d>& vRots, std::vector<Eigen::Vector3d>& vTransRot) {
  //Compute quaternions for current frame
  if(quaternions_available)
    lQuats = slerp(vIdQuat, vQuat_, (double)l/(L-1));
  else { //if quats not available, infer them
    auto bl = Matrot_.block(3*l*K, 0, 3 * K, 3); //retrieve current frame's rotation matrices
    for(int i = 0; i < K; i++) {
      Eigen::Matrix3d mat = bl.block(3*i, 0, 3, 3);
      lQuats.push_back(Eigen::Quaterniond(mat));
    }
  }

  //Compute rotation matrices for current frame
  std::vector<Eigen::Matrix3d> lRots;
  lRots.reserve(K);
  auto bl = Matrot_.block(3*l*K, 0, 3 * K, 3); //retrieve current frame's rotation matrices
  for(int i = 0; i < K; i++) {
    Eigen::Matrix3d mat = bl.block(3*i, 0, 3, 3);
    lRots.push_back(mat);
  }

  //Perform forward kinemtics both for rotation matrices and quaternions
  igl::forward_kinematics(Joints__, Links, Parents, lQuats, vTrans_, vQuats, vTrans); //quaternion-based FK
  forward_kinematics(Joints__, Links, Parents, lRots, vTrans_, vRots, vTransRot); //rotation-matrix-based FK

  //Compute the new joint positions Joints_
  Eigen::Vector3d root = Joints__.row(0).transpose();
  if(use_quaternions) {
    for(int i = 0; i < K; i++) {
      auto rotation = vQuats[i]; //A quaternion
      auto translation = vTrans[i];
      int joint_index = Links(i, 1);
      Joints_.row(joint_index) = (rotation * (Joints__.row(joint_index).transpose() - root) + translation).transpose(); //Apply affine transformation
    }
  } else {
    for(int i = 0; i < K; i++) {
      auto rotation = vRots[i]; //A rotation matrix
      auto translation = vTransRot[i];
      int joint_index = Links(i, 1);
      Joints_.row(joint_index) = (rotation * (Joints__.row(joint_index).transpose() - root) + translation).transpose(); //Apply affine transformation
    }
  }

  //Set the joint and link data
  viewer.data().set_edges(Joints_, Links, link_color);
  viewer.data().set_points(Joints_, joint_color);
}

//Update the mesh vertices, with different skinning modes
void update_mesh(vecQ& vQuats, std::vector<Eigen::Vector3d>& vTrans,
  std::vector<Eigen::Matrix3d>& vRots, std::vector<Eigen::Vector3d>& vTransRot,
  Eigen::RowVectorXd& weights) {
  //Deform mesh
  switch (skinning_mode) {
  case PER_VERTEX:
    for(int i = 0; i < V.rows(); i++) {
      Eigen::Vector3d vec = Eigen::Vector3d::Zero();
      Eigen::Vector3d vertex = V.row(i).transpose();
      if(example_data_available && context_aware) {
        if(per_vertex_unposing) {
          for(int j = 0; j < J; j++) { //for each example
            vertex += weights(j) * Displacements_[j].row(i);
          }
        }
      }
      for(int k = 0; k < K ; k++) { //for each link, compue transformation contribution
        if(use_quaternions) {
          auto rotation = vQuats[k];
          auto translation = vTrans[k];
          vec += HarmonicWeights_(i, k) * (rotation * vertex + translation).transpose();
        } else {
          auto rotation = vRots[k];
          auto translation = vTransRot[k];
          vec += HarmonicWeights_(i, k) * (rotation * vertex + translation).transpose();
        }
      }
      V_.row(i) = vec;
    }
    break;
  case DUAL_QUATERNION:
    igl::dqs(V, HarmonicWeights_, vQuats, vTrans, V_);
    break;
  case PER_FACE:
  {
    Eigen::MatrixXd Q_Tilde(3 * F.rows(), 3);
    for(int i = 0; i < F.rows(); i++) {
      Eigen::Quaterniond quat;
      blend_quaternions(vQuats, FaceWeights_.row(i), quat);
      if(example_data_available && context_aware) {
        if(!per_vertex_unposing) {
          if(equation_10) {
            //Compute the Skew and rotation matrices
            Eigen::Matrix3d def_mat_S = Eigen::Matrix3d::Zero(3, 3);
            Eigen::Matrix3d def_mat_Q = Eigen::Matrix3d::Zero(3, 3);
            vecQ quats;
            for(int j = 0; j < J; j++) { //for each example
              Eigen::Matrix3d T_jf = Deformations_[j].block(3*i, 0, 3, 3);
              Eigen::Matrix<double, 3, 3> Q, S;
              igl::polar_dec(T_jf, Q, S);
              def_mat_S += weights(j) * S;
              quats.push_back(Eigen::Quaterniond(Q));
            }
            //Transform to quaternion representation
            Eigen::Quaterniond def_quat_S = Eigen::Quaterniond(def_mat_S);
            Eigen::Quaterniond def_quat_Q;
            blend_quaternions(quats, weights, def_quat_Q);
            //Update the rotation quaternion
            quat = quat * def_quat_Q.normalized() * def_quat_S.normalized();
            quat.normalize();
          } else {
            Eigen::Matrix3d def_mat = Eigen::Matrix3d::Zero(3, 3);
            for(int j = 0; j < J; j++) {
              Eigen::Matrix3d T_jf = Deformations_[j].block(3*i, 0, 3, 3);
              def_mat += weights(j) * T_jf;
            }
            Eigen::Quaterniond def_quat = Eigen::Quaterniond(def_mat);
            quat = quat * def_quat.normalized();
            quat.normalize();
          }
        }
      }

      //Use the obtained quaternion to transform each vertex of the face:
      Eigen::Vector3d v0_prime = (quat * V.row(F(i, 0)));
      Eigen::Vector3d v1_prime = (quat * V.row(F(i, 1)));
      Eigen::Vector3d v2_prime = (quat * V.row(F(i, 2)));

      Eigen::Matrix3Xd s1(3, 3);
      s1.col(0) = (v0_prime - v2_prime).transpose();
      s1.col(1) = (v1_prime - v2_prime).transpose();
      s1.col(2) = s1.col(0).cross(s1.col(1)).normalized();

      Eigen::Matrix3Xd s2(3, 3);
      s2.col(0) = (V.row(F(i, 0)) - V.row(F(i, 2))).transpose();
      s2.col(1) = (V.row(F(i, 1)) - V.row(F(i, 2))).transpose();
      s2.col(2) = s2.col(0).cross(s2.col(1)).normalized();

      //Create the 3x3 transformation matrix
      Eigen::MatrixXd q_tilde_j = s1 * s2.inverse();
      Eigen::MatrixXd q_tilde_j_T = q_tilde_j.transpose();
      Q_Tilde.row(i) = q_tilde_j_T.row(0);
      Q_Tilde.row(i + F.rows()) = q_tilde_j_T.row(1);
      Q_Tilde.row(i + 2 * F.rows()) = q_tilde_j_T.row(2);
    }

    //Solve the Poisson-Stitching linear system and update the vertices
    Eigen::MatrixXd V_non_root = faceLBSSolver.solve(B_f * Q_Tilde - M_fc * V_root);
    igl::slice_into(V_root, root_vertices, 1, V_);
    igl::slice_into(V_non_root, non_root_vertices, 1, V_);
  }
    break;
  default:
    break;
  }
  viewer.data().set_mesh(V_, F);
}

//render each frame of the animation
void render_frame() {
  vecQ lQuats; //Relative quaternions for current frame
  vecQ vQuats; //Absolute quaternions
  std::vector<Eigen::Vector3d> vTrans; //Absolute translations obtained using base quaternions
  std::vector<Eigen::Matrix3d> vRots; //Absolute rotation matrices
  std::vector<Eigen::Vector3d> vTransRot; //Absolute translations obtained using base rotations
  Eigen::RowVectorXd weights; //example weights
  update_rig(lQuats, vQuats, vTrans, vRots, vTransRot);
  if(example_data_available && context_aware) {
    compute_example_weights(lQuats, example_poses, weights);
  }
  update_mesh(vQuats, vTrans, vRots, vTransRot, weights);
}

//Compute what to display on the next frame
bool callback_pre_draw(Viewer& viewer) {
  viewer.core().is_animating = animate;
  highlight_handles();

  if(visualize_pose) { //We visualize a given pose instead of animating
    if(example_data_available) {
      if(per_vertex_unposing)
        viewer.data().set_mesh(unposed_examples_per_vertex[current_pose], F);
      else
        viewer.data().set_mesh(unposed_examples_per_face[current_pose], F);
    } else {
      cout << "Trying to visualize unposed examples but no example is available" << endl;
    }
    return false;
  }

  if(animate) {
    //Update frame count
    if(!bidirectional) {
      ++l;
      l %= L;
    }
    else {
      if(forward_animation) {
        ++l;
        if(l == L) {
          forward_animation = false;
          l = L - 1;
        }
      } else {
        --l;
        if(l == -1) {
          forward_animation = true;
          l = 0;
        }
      }
    }
    render_frame();
  }

  return false;
}

bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers) {
  return false;
}

bool callback_mouse_down(Viewer& viewer, int button, int modifier) {
  return false;
}

bool callback_mouse_move(Viewer& viewer, int mouse_x, int mouse_y) {
  return false;
}

bool callback_mouse_up(Viewer& viewer, int button, int modifier) {
  return false;
}

int main(int argc, char *argv[]) {
  joint_color << PURPLE;
  link_color << GREENISH;
  folder = "hand";

  if(argc != 2) {
    cout << "Usage: ./assignment6 <folder_name> (i.e. <hand> or <context-aware>)" << endl;
    cout << "No folder name was specified, so using the 'hand' folder by default" << endl;
  } else {
    folder = argv[1];
    if(folder != "hand" && folder != "context-aware") {
      cout << "Supported data folders are:\n\thand\n\tcontext-aware" << endl;
      cout << "Usage: ./assignment6 <folder_name> (i.e. <hand> or <context-aware>)" << endl;
      return -1; 
    }
  }
  load_data();

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  menu.callback_draw_viewer_menu = [&]() {
    draw_viewer_menu_minimal(&viewer);

    if (ImGui::CollapsingHeader("Skeletal Animation", ImGuiTreeNodeFlags_DefaultOpen)) {
      if(ImGui::Checkbox("Animate", &animate)) {
        if(animate) visualize_pose = false;
      }
      if(ImGui::Checkbox("Use reference handles", &use_reference_handles)) {
        handle_dependent_processing();
      }
      if(ImGui::SliderInt("Frame", &l, 0, L-1, "%d", 0)) {
        render_frame();
      }
      ImGui::Checkbox("Quaternion-based", &use_quaternions);
      ImGui::Checkbox("Bidirectional", &bidirectional);
    }

    if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::Combo("Handle color mode", (int*)(&handles_color_mode), "All handles\0Harmonic Single Handle\0Single Handle\0\0");
      if(handles_color_mode != ALL_HANDLES) {
        ImGui::SliderInt("Joint index", &selected_handle, 0, K-1, "%d", 0);
      }
      if(ImGui::Combo("Skinning Mode", (int*)(&skinning_mode), "Per-Vertex Linear\0Dual Quaternion\0Per-Face Linear\0\0")) {
        per_vertex_unposing = skinning_mode != PER_FACE;
      }
    }

    if (ImGui::CollapsingHeader("Context-based", ImGuiTreeNodeFlags_DefaultOpen)) {
      ImGui::Checkbox("Context-aware", &context_aware);
      if(context_aware) {
        if(!per_vertex_unposing) {
          ImGui::Checkbox("Use Equation 10", &equation_10);
        }
        ImGui::Combo("RBF Kernel", (int*)(&rbf_function), "Gaussian\0Polyharmonic r^3\0Thin plate\0\0");
        if(rbf_function == GAUSSIAN) {
          ImGui::InputDouble("Epsilon", &epsilon);
        }
      }
      if(ImGui::Checkbox("Visualize pose", &visualize_pose)) {
        if(visualize_pose) animate = false;
      }
      if(visualize_pose) {
        ImGui::SliderInt("Pose index", &current_pose, 0, J-1, "%d", 0);
      }
      if(ImGui::Checkbox("Per-vertex unposing", &per_vertex_unposing)) {
        if(per_vertex_unposing)
          skinning_mode = PER_VERTEX;
        else
          skinning_mode = PER_FACE;
      }
    }
  };

  viewer.callback_key_down = callback_key_down;
  viewer.callback_mouse_down = callback_mouse_down;
  viewer.callback_mouse_move = callback_mouse_move;
  viewer.callback_mouse_up = callback_mouse_up;
  viewer.callback_pre_draw = callback_pre_draw;

  viewer.data().point_size = 5;
  viewer.data().show_overlay_depth = false;
  viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);
  viewer.launch();
}

//Reduced version of the ImGUI::draw_viewer_menu
void draw_viewer_menu_minimal(Viewer* viewer)
{
  // Viewing options
  if (ImGui::CollapsingHeader("Viewing Options", ImGuiTreeNodeFlags_DefaultOpen))
  {
    if (ImGui::Button("Snap canonical view", ImVec2(-1, 0)))
    {
      viewer->snap_to_canonical_quaternion();
    }

    ImGui::PushItemWidth(80);

    // Select rotation type
    int rotation_type = static_cast<int>(viewer->core().rotation_type);
    static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
    static bool orthographic = true;
    if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0"))
    {
      using RT = igl::opengl::ViewerCore::RotationType;
      auto new_type = static_cast<RT>(rotation_type);
      if (new_type != viewer->core().rotation_type)
      {
        if (new_type == RT::ROTATION_TYPE_NO_ROTATION)
        {
          trackball_angle = viewer->core().trackball_angle;
          orthographic = viewer->core().orthographic;
          viewer->core().trackball_angle = Eigen::Quaternionf::Identity();
          viewer->core().orthographic = true;
        }
        else if (viewer->core().rotation_type == RT::ROTATION_TYPE_NO_ROTATION)
        {
          viewer->core().trackball_angle = trackball_angle;
          viewer->core().orthographic = orthographic;
        }
        viewer->core().set_rotation_type(new_type);
      }
    }

    // Orthographic view
    ImGui::Checkbox("Orthographic view", &(viewer->core().orthographic));
    ImGui::PopItemWidth();
  }

  // Helper for setting viewport specific mesh options
  auto make_checkbox = [&](const char *label, unsigned int &option)
  {
    return ImGui::Checkbox(label,
      [&]() { return viewer->core().is_set(option); },
      [&](bool value) { return viewer->core().set(option, value); }
    );
  };

  // Draw options
  if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen))
  {
    make_checkbox("Show overlay", viewer->data().show_overlay);
  }

  // Overlays
  if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen))
  {
    make_checkbox("Wireframe", viewer->data().show_lines);
    make_checkbox("Fill", viewer->data().show_faces);
  }
}