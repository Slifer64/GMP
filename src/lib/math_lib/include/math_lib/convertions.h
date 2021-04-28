#ifndef MATH_LIB_CONVERTIONS_H
#define MATH_LIB_CONVERTIONS_H

#include <ros/ros.h>
#include <geometry_msgs/Transform.h>
#include <math_lib/math.h>

namespace as64_
{

namespace math_
{

void rosTransform_to_eigenTransform(const geometry_msgs::Transform &rosTrans, Eigen::Matrix4d &eigenTrans);

} // namespace math_

} // namespace as64_

#endif // MATH_LIB_CONVERTIONS_H
