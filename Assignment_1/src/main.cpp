////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
// Shortcut to avoid  everywhere, DO NOT USE IN .h
using namespace Eigen;
////////////////////////////////////////////////////////////////////////////////

const std::string root_path = DATA_DIR;

// Computes the determinant of the matrix whose columns are the vector u and v
double inline det(const Vector2d &u, const Vector2d &v)
{
    return u.x() * v.y() - u.y() * v.x();
}

// Return true iff [a,b] intersects [c,d]
bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{
    // TODO
    return (det(b - a, c - a) * det(b - a, d - a) < 0 &&
            det(d - c, a - c) * det(d - c, b - c) < 0);
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    // TODO
    Vector2d outside(std::numeric_limits<double>::min(), std::numeric_limits<double>::min());
    for (auto point : poly)
    {
        outside.x() = std::max(outside.x(), point.x());
        outside.y() = std::max(outside.y(), point.y());
    }
    outside.x() += 1; outside.y() -= 1;
    // 2. Cast a ray from the query point to the 'outside' point, count number of intersections
    int intersection = 0;
    for (int i = 0; i < poly.size(); ++i)
    {
        if (intersect_segment(query, outside, poly[i], poly[(i + 1) % poly.size()]))
        {
            intersection++;
        }
    }
    return intersection % 2 == 1;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Vector2d> load_xyz(const std::string &filename)
{
    std::vector<Vector2d> points;
    std::ifstream in(filename);

    double x; double y; double z;
    std::string line;

    while (getline(in, line)) {
        in >> x >> y >> z;
        points.push_back(Vector2d(x, y));
    }

    in.close();
    return points;
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
    std::ofstream out(filename);
    out << points.size() << std::endl;
    for (size_t i = 0; i < points.size(); ++i)
    {
        out << points[i].x() << " " << points[i].y() << " 0" << std::endl;
    }
    out.close();
}

std::vector<Vector2d> load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    std::vector<Vector2d> points;
    std::vector<Vector2d> poly;
    char key;
    while (in >> key)
    {
        if (key == 'v')
        {
            double x, y, z;
            in >> x >> y >> z;
            points.push_back(Vector2d(x, y));
        }
        else if (key == 'f')
        {
            std::string line;
            std::getline(in, line);
            std::istringstream ss(line);
            int id;
            while (ss >> id)
            {
                poly.push_back(points[id - 1]);
            }
        }
    }
    return poly;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Vector2d> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    std::vector<Vector2d> poly = load_obj(poly_path);
    std::vector<Vector2d> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);

    return 0;
}
