// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"


// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;



void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

Vector3d ray_parallelogram_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction,
                                    const Vector3d &origin, const Vector3d &u, const Vector3d &v)
{

    Matrix3d m_matrix;
    Vector3d b_vector = ray_origin - origin;
    m_matrix << u, v, -ray_direction;
    
    return m_matrix.inverse() * b_vector;
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            Vector3d intersection = ray_parallelogram_intersection(ray_origin, ray_direction,
                                                                   pgram_origin, pgram_u, pgram_v);
            double u = intersection(0);
            double v = intersection(1);
            double t = intersection(2);

            // Check if the ray intersects with the parallelogram
            if (t > 0 && u >= 0 && u <= 1 && v >= 0 && v <= 1)
            {

                Vector3d ray_intersection = ray_origin + t * ray_direction;

                Vector3d ray_normal = -(pgram_u.cross(pgram_v).normalized());

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}


void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray (origin point and direction)
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = pixel_center - camera_origin;

            Vector3d intersection = ray_parallelogram_intersection(ray_origin, ray_direction,
                                                                   pgram_origin, pgram_u, pgram_v);
            double u = intersection(0);
            double v = intersection(1);
            double t = intersection(2);

            // Check if the ray intersects with the parallelogram
            if (t > 0 && u >= 0 && u <= 1 && v >= 0 && v <= 1)
            {
                // The ray hit the parallelogram, compute the exact intersection point
                Vector3d ray_intersection = ray_origin + t * ray_direction;

                // Compute normal at the intersection point
                Vector3d ray_normal = -(pgram_u.cross(pgram_v).normalized());

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");

    MatrixXd R = MatrixXd::Zero(800, 800);
    MatrixXd G = MatrixXd::Zero(800, 800);
    MatrixXd B = MatrixXd::Zero(800, 800);


    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray (origin point and direction)
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = pixel_center - camera_origin;


            double h = pow((ray_direction.dot(ray_origin - sphere_center)), 2) -
                          ((ray_direction.dot(ray_direction)) *
                           ((ray_origin - sphere_center).dot(ray_origin - sphere_center) - pow(sphere_radius, 2)));

            // Intersect with the sphere
            // implement the generic ray sphere intersection
            if (h >= 0)
            {
                double t_distance = (-ray_direction.dot(ray_origin - sphere_center) - sqrt(h)) / ray_direction.dot(ray_direction);
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection = ray_origin + t_distance * ray_direction;

                // Compute normal at the intersection point
                Vector3d ray_normal = (2*(ray_intersection - sphere_center)).normalized();


                Vector3d h = (light_position - ray_intersection).normalized() + (camera_origin - ray_intersection).normalized();

                // Add shading parameter here
                const double diffuse = std::max(0., (light_position - ray_intersection).normalized().dot(ray_normal));
                const double specular = std::max(0., h.normalized().dot(ray_normal));

                // Simple diffuse model
                R(i, j) = ambient + diffuse_color(0) * diffuse + specular_color(0) * pow(specular, specular_exponent);
                G(i,j) = ambient + diffuse_color(1) * diffuse + specular_color(1) * pow(specular, specular_exponent);
                B(i,j) = ambient + diffuse_color(2) * diffuse + specular_color(2) * pow(specular, specular_exponent);


                // Clamp to zero
                R(i, j) = std::max(R(i, j), 0.);
                G(i, j) = std::max(G(i, j), 0.);
                B(i, j) = std::max(B(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

int main()
{
//    raytrace_sphere();
//    raytrace_parallelogram();
//    raytrace_perspective();
    raytrace_shading();

    return 0;
}
