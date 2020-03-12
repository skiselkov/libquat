/*
 * CONFIDENTIAL
 *
 * Copyright 2020 Saso Kiselkov. All rights reserved.
 *
 * NOTICE:  All information contained herein is, and remains the property
 * of Saso Kiselkov. The intellectual and technical concepts contained
 * herein are proprietary to Saso Kiselkov and may be covered by U.S. and
 * Foreign Patents, patents in process, and are protected by trade secret
 * or copyright law. Dissemination of this information or reproduction of
 * this material is strictly forbidden unless prior written permission is
 * obtained from Saso Kiselkov.
 */

#include <stdio.h>

#include <acfutils/assert.h>
#include <acfutils/perf.h>

#include "libquat.h"

/*
 * Calculates the Hamilton product of quaternions `p' and `q':
 * H(p,q) = q.p.q^-1
 */
struct quat
quat_hamil_prod(struct quat p, struct quat q)
{
	struct quat q_inv, out;

	ASSERT(!IS_NULL_QUAT(p));
	ASSERT(!IS_NULL_QUAT(q));
	quat_inverse(q_inv.v, q.v);
	quat_multiply(out.v, q.v, p.v);
	quat_multiply(out.v, out.v, q_inv.v);

	return (out);
}

/*
 * Constructs a quaternion to translate from X-Plane OpenGL local coordinates
 * with origin at `refpt' to the Earth-Centered-modified-Inertial frame at
 * `ref_time'.
 */
struct quat
quat_local2ecmi(geo_pos2_t refpt, double ref_time)
{
	struct quat lat_q, lon_q, out_q;

	ASSERT(isfinite(ref_time));
	ASSERT(!IS_NULL_GEO_POS(refpt));
	quat_from_axis_angle(lat_q.v, QUAT_AXIS_X.v, DEG2RAD(90 - refpt.lat));
	quat_from_axis_angle(lon_q.v, QUAT_AXIS_Y.v,
	    DEG2RAD(refpt.lon + EARTH_ROT_RATE * ref_time));
	quat_multiply(out_q.v, lon_q.v, lat_q.v);

	return (out_q);
}

/*
 * Constructs inverse quaternion to translate from ECmI to X-Plane
 * OpenGL local coordinates.
 */
struct quat
quat_ecmi2local(geo_pos2_t refpt, double ref_time)
{
	struct quat q, q_inv;

	ASSERT(isfinite(ref_time));
	ASSERT(!IS_NULL_GEO_POS(refpt));
	q = quat_local2ecmi(refpt, ref_time);
	quat_inverse(q_inv.v, q.v);

	return (q_inv);
}

/*
 * Constructs a relative rotation quaternion that expresses the relative
 * rotation from `q1' to `q2'. This can subsequently be concatenated onto
 * `q1' to obtain `q2'.
 */
struct quat
quat_rot_rel(struct quat q1, struct quat q2)
{
	struct quat q1_inv, d_quat;

	ASSERT(!IS_NULL_QUAT(q1));
	ASSERT(!IS_NULL_QUAT(q2));
	quat_inverse(q1_inv.v, q1.v);
	quat_multiply(d_quat.v, q2.v, q1_inv.v);

	return (d_quat);
}

/*
 * Concatenates a relative rotation quaternion `delta' onto `from'.
 * This is equivalent to quaternion multiplication of `delta'.`from'.
 */
struct quat
quat_rot_concat(struct quat from, struct quat delta)
{
	struct quat to;

	ASSERT(!IS_NULL_QUAT(from));
	ASSERT(!IS_NULL_QUAT(delta));
	quat_multiply(to.v, delta.v, from.v);

	return (to);
}

void
quat_to_axis_angle(struct quat q, struct vec3 *axis, double *angle)
{
	double half = acos(q.w);
	double s = sin(half);

	ASSERT(!IS_NULL_QUAT(q));
	/*
	 * Zero roration quaternion has an arbitrary axis, so just pick one.
	 */
	if (s == 0) {
		*axis = QUAT_AXIS_X;
		*angle = 0;
		return;
	}
	axis->x = q.x / s;
	axis->y = q.y / s;
	axis->z = q.z / s;
	*angle = 2 * half;
}

struct quat
quat_from_euler(double psi, double theta, double phi)
{
	struct quat q;
	double sin_psi = sin(psi / 2), cos_psi = cos(psi / 2);
	double sin_theta = sin(theta / 2), cos_theta = cos(theta / 2);
	double sin_phi = sin(phi / 2), cos_phi = cos(phi / 2);

	ASSERT(isfinite(psi));
	ASSERT(isfinite(theta));
	ASSERT(isfinite(phi));
	q.w = cos_phi * cos_theta * cos_psi + sin_phi * sin_theta * sin_psi;
	q.x = sin_phi * cos_theta * cos_psi - cos_phi * sin_theta * sin_psi;
	q.y = cos_phi * sin_theta * cos_psi + sin_phi * cos_theta * sin_psi;
	q.z = cos_phi * cos_theta * sin_psi - sin_phi * sin_theta * cos_psi;

	return (q);
}

void
quat_to_euler(struct quat q, double *psi, double *theta, double *phi)
{
	double sq_w = q.w * q.w;
	double sq_x = q.x * q.x;
	double sq_y = q.y * q.y;
	double sq_z = q.z * q.z;

	ASSERT(!IS_NULL_QUAT(q));
	*psi = atan2(2.0 * (q.x * q.y + q.z * q.w),
	    (sq_x - sq_y - sq_z + sq_w));
	*theta = asin(-2.0 * (q.x * q.z - q.y * q.w));
	*phi = atan2(2.0 * (q.y * q.z + q.x * q.w),
	    (-sq_x - sq_y + sq_z + sq_w));
}

/*
 * Converts a vect3_t in OpenGL space to a quaternion suitable for
 * transformation using Hamilton products.
 */
struct quat
quat_from_vect3_gl(vect3_t v)
{
	struct quat q;
	quat(q.v, -v.z, v.x, -v.y, 0);
	return (q);
}

/*
 * Converts a quat back to a vect3_t in OpenGL space.
 */
vect3_t
quat_to_vect3_gl(struct quat q)
{
	return (VECT3(q.x, -q.z, -q.y));
}

void
quat_print(struct quat q)
{
	printf("%11f %11f %11f %11f\n", q.w, q.x, q.y, q.z);
}

void
quat_print_euler(struct quat q)
{
	double psi, theta, phi;

	ASSERT(!IS_NULL_QUAT(q));
	quat_to_euler(q, &psi, &theta, &phi);
	printf("psi: %f  theta: %f  phi: %f\n",
	    RAD2DEG(psi), RAD2DEG(theta), RAD2DEG(phi));
}

void
quat_print_axis_angle(struct quat q)
{
	struct vec3 axis;
	double angle;

	ASSERT(!IS_NULL_QUAT(q));
	quat_to_axis_angle(q, &axis, &angle);
	printf("%11f %11f %11f %11f deg\n", axis.y, -axis.z, -axis.x,
	    RAD2DEG(angle));
}