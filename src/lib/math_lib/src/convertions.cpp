#include <math_lib/convertions.h>

namespace as64_
{

namespace math_
{

void rosTransform_to_eigenTransform(const geometry_msgs::Transform &rosTrans, Eigen::Matrix4d &eigenTrans)
{
	eigenTrans(0,3) = rosTrans.translation.x;
	eigenTrans(1,3) = rosTrans.translation.y;
	eigenTrans(2,3) = rosTrans.translation.z;

	Eigen::Vector4d quat;
	quat << rosTrans.rotation.w, rosTrans.rotation.x, rosTrans.rotation.y, rosTrans.rotation.z;

	eigenTrans.block(0,0,3,3) = as64_::math_::quat2rotm(quat);

	eigenTrans.row(3) << 0, 0, 0, 1;
}

} // namespace math_

} // namespace as64_
