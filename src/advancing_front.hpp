#pragma once
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include <nanoflann.hpp>

template <int dim>
using Point = Eigen::Matrix<double, dim, 1>;

template <int dim>
using Normal = Eigen::Matrix<double, dim, 1>;

template <int dim>
struct Nodes {
    std::vector<Point<dim>> points;
    std::vector<Normal<dim>> normals;
    
    Nodes() = default;

    Nodes(size_t Npts, double *p) : points(Npts)
    {
        for (size_t i = 0; i < Npts; i++) {
            for (int j = 0; j < dim; j++) {
                points[i][j] = p[i + j*Npts];
            }
        }
    }
    
    Nodes(size_t Npts, double *p, double *n) : points(Npts), normals(Npts)
    {
        for (size_t i = 0; i < Npts; i++) {
            for (int j = 0; j < dim; j++) {
                points[i][j] = p[i + j*Npts];
                normals[i][j] = n[i + j*Npts];
            }
        }
    }

    inline size_t kdtree_get_point_count() const
    {
        return points.size();
    }

    inline double kdtree_get_pt(const size_t i, int j) const
    {
        return points[i][j];
    }
    
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const
    {
        return false;
    }
};

template <int dim>
using StaticKDTree = nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, Nodes<dim>>, Nodes<dim>, dim, size_t>;

template <int dim>
using DynamicKDTree = nanoflann::KDTreeSingleIndexDynamicAdaptor<
        nanoflann::L2_Simple_Adaptor<double, Nodes<dim>>, Nodes<dim>, dim, size_t>;

template <int dim>
using AABB = Eigen::AlignedBox<double,dim>;

template <int dim>
Point<dim> aabb_random_point(const AABB<dim> &aabb, std::mt19937_64 &rng)
{
    Point<dim> p;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    for (size_t i = 0; i < dim; i++) {
        p[i] = uniform(rng);
    }
    return aabb.min() + p.cwiseProduct(aabb.diagonal());
}

template <int dim>
Eigen::Matrix<double,dim,dim> get_random_rotation(std::mt19937_64 &rng)
{
    static_assert(1 <= dim && dim <= 3, "Invalid dimension");
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    Eigen::Matrix<double,dim,dim> random_rotation;
    if constexpr (dim == 1) {
        random_rotation(0,0) = 1.0;
    } else if constexpr (dim == 2) {
        double theta = 2*M_PI*uniform(rng);
        random_rotation(0,0) = cos(theta);
        random_rotation(1,0) = sin(theta);
        random_rotation(0,1) = -sin(theta);
        random_rotation(1,1) = cos(theta);
    } else if constexpr (dim == 3) {
        double a = uniform(rng);
        double b = uniform(rng);
        double c = uniform(rng);
        double w = sqrt(1-a)*sin(2*M_PI*b);
        double x = sqrt(1-a)*cos(2*M_PI*b);
        double y =   sqrt(a)*sin(2*M_PI*c);
        double z =   sqrt(a)*cos(2*M_PI*c);
        Eigen::Quaternion<double> q(w, x, y, z);
        random_rotation = q.toRotationMatrix();
    }
    return random_rotation;
}

double candidates_buf_1D[] = {
    1.0,
   -1.0
};

double candidates_buf_2D[] = {
    1.000000000000000e+00,                         0,
    9.135454576426009e-01,     4.067366430758002e-01,
    6.691306063588582e-01,     7.431448254773942e-01,
    3.090169943749475e-01,     9.510565162951535e-01,
   -1.045284632676533e-01,     9.945218953682734e-01,
   -4.999999999999998e-01,     8.660254037844387e-01,
   -8.090169943749473e-01,     5.877852522924732e-01,
   -9.781476007338057e-01,     2.079116908177593e-01,
   -9.781476007338057e-01,    -2.079116908177591e-01,
   -8.090169943749475e-01,    -5.877852522924730e-01,
   -5.000000000000004e-01,    -8.660254037844385e-01,
   -1.045284632676542e-01,    -9.945218953682733e-01,
    3.090169943749472e-01,    -9.510565162951536e-01,
    6.691306063588585e-01,    -7.431448254773940e-01,
    9.135454576426010e-01,    -4.067366430758002e-01
};

double candidates_buf_3D[] = {
    2.769940904098851e-01,     3.097530659272069e-01,    -9.095753470860432e-01,
   -9.735855355968755e-02,     5.138546962740651e-01,    -8.523348304310071e-01,
    3.263744552205621e-01,     6.490754269993295e-01,    -6.871541348526685e-01,
    6.566366526278439e-01,     1.712256625006458e-01,    -7.345134981243812e-01,
    1.667126870118809e-01,     1.836084708895852e-02,    -9.858345496499169e-01,
   -2.499679908546314e-01,     2.238943686622654e-01,    -9.420123753058791e-01,
   -4.816717831919198e-01,     5.753002487947900e-01,    -6.610763322138896e-01,
   -8.142239235418715e-02,     8.080285986136718e-01,    -5.834896553031137e-01,
    2.589787495372250e-01,     8.916801141717987e-01,    -3.712634930594486e-01,
    6.480154223000650e-01,     6.324355178684952e-01,    -4.243834683392799e-01,
    8.917455417736815e-01,     1.833987276392105e-01,    -4.137085875674832e-01,
   -1.366803346519356e-01,    -4.295451264830512e-02,    -9.896834827168641e-01,
   -6.304319980036653e-01,     1.593967891696202e-01,    -7.597026783522104e-01,
   -7.958071393494702e-01,     4.300386623859934e-01,    -4.263305593241921e-01,
   -5.064825784761433e-01,     7.972912125186770e-01,    -3.283323318539545e-01,
   -7.280144976556017e-02,     9.852161194727271e-01,    -1.550778734801765e-01,
    3.721699847134812e-01,     9.281210458176063e-01,     8.990372005643484e-03,
    7.261176686104671e-01,     6.864900460375882e-01,    -3.852983289620796e-02,
    9.622254696252533e-01,     2.703912518251204e-01,    -3.179176844569755e-02,
   -2.917007022073100e-01,    -4.156522215545622e-01,    -8.614777600428924e-01,
   -7.392636383657355e-01,    -2.312224658080010e-01,    -6.324756471959397e-01,
   -8.241830645064856e-01,     5.639520137661819e-01,    -5.177260230823361e-02,
   -9.303247863123811e-01,     5.584803063862490e-02,    -3.624593624761395e-01,
   -4.924413134881133e-01,     8.680892239123101e-01,     6.263107932508237e-02,
   -1.182874449496023e-02,     9.713514404477581e-01,     2.373530280905182e-01,
    3.841287841366001e-01,     8.346101196560495e-01,     3.948050472897376e-01,
    7.373937251932121e-01,     5.812547823173372e-01,     3.440833795447008e-01,
    9.184775940250842e-01,     1.797462835527348e-01,     3.522643649631231e-01,
    8.036168845943922e-02,    -4.641972383186805e-01,    -8.820787509997382e-01,
   -1.955142546812522e-01,    -7.514771112418902e-01,    -6.301240572268115e-01,
   -5.938577747017016e-01,    -6.040658247174705e-01,    -5.314484197311581e-01,
   -9.717137709820451e-01,     2.207253452108027e-01,     8.398017185875187e-02,
   -7.689627124460980e-01,     5.407428054280582e-01,     3.410184236156938e-01,
   -9.628318642800963e-01,    -2.512101861155437e-01,    -9.923831678695647e-02,
   -4.362815472018196e-01,     7.813344926537603e-01,     4.462900650481450e-01,
    2.875701477957233e-02,     8.034179086844044e-01,     5.947206891526016e-01,
    3.627964225920945e-01,     6.010896270843681e-01,     7.120884888593226e-01,
    6.988997527607246e-01,     2.482247754113347e-01,     6.707634429983425e-01,
    4.985392233092584e-01,    -4.601508066829655e-01,    -7.346562991843311e-01,
    1.739438094631374e-01,    -8.643220696677967e-01,    -4.719013785046916e-01,
   -2.726714089532511e-01,    -9.208237776357632e-01,    -2.788079505323539e-01,
   -7.180154018637394e-01,    -6.771767626677451e-01,    -1.608897597402811e-01,
   -8.949518251509300e-01,     6.893640486904842e-02,     4.408049486368660e-01,
   -6.614761821190372e-01,     3.454484581309826e-01,     6.656685536091885e-01,
   -8.631180230733720e-01,    -4.495903115223914e-01,     2.299909346715935e-01,
   -3.355077064885693e-01,     5.581547311522115e-01,     7.588793546929462e-01,
    2.190771009792657e-01,     3.024007449538551e-01,     9.276632003469115e-01,
    6.506775727282067e-01,    -6.372680792003208e-01,    -4.129262544097470e-01,
    4.743283296519485e-01,    -1.286908121681012e-01,    -8.708911014318076e-01,
    3.453274037989134e-01,    -9.313190817324873e-01,    -1.157313794372931e-01,
   -5.553644809845049e-02,    -9.972305713342915e-01,     4.946605430889624e-02,
   -4.763915062658228e-01,    -8.747384052574062e-01,     8.879107570871300e-02,
   -7.651338166059831e-01,    -2.567639313641796e-01,     5.904595889951934e-01,
   -5.773016731654850e-01,     2.276913668490853e-02,     8.162134185217463e-01,
   -2.599258996588312e-01,     2.351162193196809e-01,     9.365676110667969e-01,
    6.412789356824490e-02,     1.888329081033493e-02,     9.977630152469433e-01,
    8.551627261361200e-01,    -2.986478216959922e-01,    -4.236817088613556e-01,
    8.013577508829152e-01,    -5.964866726614431e-01,    -4.504891160900738e-02,
    5.434825039156430e-01,    -8.132452508459997e-01,     2.079878119362043e-01,
    7.156882710638023e-02,    -9.076770864422293e-01,     4.135217161581273e-01,
   -3.529955391834372e-01,    -8.143265265997235e-01,     4.607238407035436e-01,
   -5.032596694676454e-01,    -5.016583978861989e-01,     7.036110835664603e-01,
   -3.446689159652604e-01,    -2.181581377483625e-01,     9.130226532247153e-01,
    9.976320828131545e-02,    -2.875538683602356e-01,     9.525544997870112e-01,
    9.449806157625878e-01,    -2.869150147221864e-01,     1.571350061568965e-01,
    6.601381520709284e-01,    -5.703530124695579e-01,     4.887893834232491e-01,
    2.577149864229281e-01,    -6.942253800293761e-01,     6.720372813290186e-01,
   -1.324871834979424e-01,    -6.600134597556269e-01,     7.394791269198815e-01,
    8.285407688008155e-01,    -1.945901084400290e-01,     5.250284603069163e-01,
    4.833022137344157e-01,    -3.800011969194029e-01,     7.886812160431073e-01,
    5.127840626789326e-01,    -2.424066905407739e-02,     8.581753288380523e-01
};

constexpr size_t Ncandidates[4] = {0, sizeof(candidates_buf_1D)/sizeof(double),
    sizeof(candidates_buf_2D)/(2*sizeof(double)), sizeof(candidates_buf_3D)/(3*sizeof(double))};

template <int dim>
bool inclusion_query_strict(const Point<dim> &y, const Nodes<dim> &Z, const StaticKDTree<dim> &kdtree)
{
    // Get indices of the dim nearest neighbors
    size_t Knn_index[dim];
    double Knn_dist2[dim];
    size_t K = kdtree.knnSearch(y.data(), dim, Knn_index, Knn_dist2);
    
    // Check if y is in the half-space dot(z_i-y,nz_i) > 0 for all K nearest-neighbors z_i
    for (size_t i = 0; i < K; i++) {
        auto z = Z.points[Knn_index[i]];
        auto nz = Z.normals[Knn_index[i]];
        if (nz.dot(z - y) <= 0) {
            return false;
        }
    }
    return true;
}

template <int dim>
bool inclusion_query(const Point<dim> &y, const Nodes<dim> &Z, const StaticKDTree<dim> &kdtree)
{
    // Get indices of the dim nearest neighbors
    size_t Knn_index[dim];
    double Knn_dist2[dim];
    size_t K = kdtree.knnSearch(y.data(), dim, Knn_index, Knn_dist2);
    
    // Check if y is inside or outside for all K nearest-neighbors z_i
    bool all_inside = true;
    bool any_inside = false;
    for (size_t i = 0; i < K; i++) {
        auto z = Z.points[Knn_index[i]];
        auto nz = Z.normals[Knn_index[i]];
        bool is_inside = (nz.dot(z - y) > 0);
        all_inside = all_inside && is_inside;
        any_inside = any_inside || is_inside;
    }
    if (all_inside) {
        return true;
    }
    if (!any_inside) {
        return false;
    }
    
    // Edge case: y is inside wrt to some z_i, and outside wrt to some z_j.
    // Start by rejecting points too close to the boundary.
    for (size_t i = 0; i < K; i++) {
        auto z_i = Z.points[Knn_index[i]];
        double maxr2 = 0.0;
        for (size_t j = 0; j < K; j++) {
            auto z_j = Z.points[Knn_index[j]];
            maxr2 = std::max(maxr2,(z_j-z_i).squaredNorm());
        }
        if ((y-z_i).squaredNorm() <= maxr2) {
            return false;
        }
    }
    
    // Wiggle the point y to solve the ambiguity. Include points far enough
    // in the interior, and reject points far enough in the exterior.
    double min_distance = sqrt(Knn_dist2[0]);
    for (size_t i = 0; i < dim; i++) {
        Point<dim> dy = Point<dim>::Zero();
        dy(i) = 0.5 * min_distance;
        if (!inclusion_query_strict<dim>(y + dy, Z, kdtree) ||
            !inclusion_query_strict<dim>(y - dy, Z, kdtree)) {
            return false;
        }
    }
    return true;
}

template <int d, int n>
Point<n> normal_from_jacobian(const Eigen::Matrix<double,n,d> &J)
{
    static_assert(1 <= d && d == n-1 && n <= 3, "Invalid dimensions");
    if constexpr (d == 1 && n == 2) {
        Eigen::Matrix<double,2,1> nu = {J(1), -J(0)};
        return nu.normalized();
    } else if constexpr (d == 2 && n == 3) {
        auto nu = J.col(0).cross(J.col(1));
        return nu.normalized();
    }
}

template <int d, int n, typename G_t, typename dG_t, typename h_t>
Nodes<n> advancing_front(const AABB<d> &aabb, const Nodes<d> &Z, Nodes<d> &Y,
    G_t G, dG_t dG, h_t h, size_t Nmax, std::mt19937_64 rng)
{
    static_assert(1 <= d && d <= n && n <= 3, "Invalid dimensions");
    
    // Initialize the set Y, if empty
    if (Y.points.empty()) {
        if (Z.points.empty()) {
            Y.points.push_back(aabb_random_point(aabb,rng));
        } else {
            Y.points = Z.points;
        }
    }

    // Initialize the set GY as G(Y)
    size_t NY = Y.points.size();
    Nodes<n> GY;
    for (size_t i = 0; i < NY; i++) {
        auto y = Y.points[i];
        auto Gy = G(y);
        GY.points.push_back(Gy);
        if constexpr (d == n-1) {
            auto nGy = normal_from_jacobian<d,n>(dG(y));
            GY.normals.push_back(nGy);
        }
    }
    
    // Initialize kd-trees to speed up nearest-neighbor queries
    StaticKDTree<d> kdtree_Z(d, Z, {20});
    DynamicKDTree<n> kdtree_GY(n, GY, {20});
    
    // Initialize sets of candidates
    double *candidates_buf = nullptr;
    if constexpr (d == 1) {
        candidates_buf = candidates_buf_1D;
    } else if constexpr (d == 2) {
        candidates_buf = candidates_buf_2D;
    } else {
        candidates_buf = candidates_buf_3D;
    }
    Eigen::Map<Eigen::Matrix<double,d,Ncandidates[d]>> candidates(candidates_buf);
    
    // Iterate over Y. Every point is a center for expansion.
    for (size_t i = 0; i < NY && NY < Nmax; i++) {
        // Get current center of expansion and a random rotation matrix
        const Point<d> y = Y.points[i];
        const Point<n> Gy = G(y);
        const double hGy = h(Gy);
        const Eigen::Matrix<double,n,d> dGy = dG(y);
        const Eigen::Matrix<double,d,d> random_rotation = get_random_rotation<d>(rng);
        
        // Iterate over candidate nodes around y
        for (size_t j = 0; j < Ncandidates[d]; j++) {
            // Generate new candidate x around y
            auto dy = random_rotation * candidates.col(j);
            double alpha = hGy / ((dGy*dy).norm()+1e-15);
            auto x = y + alpha * dy;
            
            // Test if x is inside the parametric domain
            if (!aabb.contains(x)) {
                continue;
            }
            if (!Z.points.empty() && !inclusion_query<d>(x,Z,kdtree_Z)) {
                continue;
            }
            
            // Test if G(x) is too close to an already existing node
            Point<n> Gx = G(x);
            double distance2_Gx_Gy = (Gx-Gy).squaredNorm();
            nanoflann::KNNResultSet<double> result_set(1);
            size_t nn_index;
            double distance2_Gx_nn;
            result_set.init(&nn_index,&distance2_Gx_nn);
            bool knn_found = kdtree_GY.findNeighbors(result_set, Gx.data());
            double reltol = 1-1e-10;
            if (knn_found && distance2_Gx_nn < distance2_Gx_Gy * reltol) {
                continue;
            }
            
            // All tests passed, we can add x to Y and G(x) to GY
            Y.points.push_back(x);
            GY.points.push_back(Gx);
            if constexpr (d == n-1) {
                auto nGx = normal_from_jacobian<d,n>(dG(x));
                GY.normals.push_back(nGx);
            }
            kdtree_GY.addPoints(NY,NY);
            NY++;
        }
    }
    return GY;
}
