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

#ifndef	_LIBQUAT_H_
#define	_LIBQUAT_H_

#include <acfutils/geom.h>
#include <mathc.h>

#ifdef __cplusplus
extern "C" {
#endif

#define	QUAT_VEC3(_x, _y, _z) ((struct vec3){ .x = (_x), .y = (_y), .z = (_z)})
#define	QUAT_AXIS_X		QUAT_VEC3(0, 1, 0)
#define	QUAT_AXIS_Y		QUAT_VEC3(0, 0, -1)
#define	QUAT_AXIS_Z		QUAT_VEC3(-1, 0, 0)

#define	IS_NULL_QUAT(q)		(!isfinite((q).w))
#define	NULL_QUAT		((struct quat){.v = {NAN, NAN, NAN, NAN}})

struct quat quat_hamil_prod(struct quat p, struct quat q);
struct quat quat_local2ecmigl(geo_pos2_t refpt, double ref_time);
struct quat quat_ecmigl2local(geo_pos2_t refpt, double ref_time);
struct quat quat_rot_rel(struct quat q1, struct quat q2);
struct quat quat_rot_concat(struct quat from, struct quat delta);
void quat_to_axis_angle(struct quat q, struct vec3 *axis, double *angle);
struct quat quat_from_euler(double psi, double theta, double phi);
struct quat quat_from_euler_deg(double psi, double theta, double phi);
void quat_to_euler(struct quat q, double *psi, double *theta, double *phi);
void quat_to_euler_deg(struct quat q, double *psi, double *theta, double *phi);
struct quat quat_from_xp(struct quat xp_q);
struct quat quat_to_xp(struct quat q);

/*
 * Converts a vect3_t in OpenGL space to a quaternion suitable for
 * transformation using Hamilton products.
 */
static inline struct quat
quat_from_vect3_gl(vect3_t v)
{
	return ((struct quat){ .x = -v.z, .y = v.x, .z = -v.y });
}

/*
 * Converts a quat back to a vect3_t in OpenGL space.
 */
static inline vect3_t
quat_to_vect3_gl(struct quat q)
{
	return (VECT3(q.y, -q.z, -q.x));
}

void quat_print(struct quat q);
void quat_print_euler(struct quat q);
void quat_print_axis_angle(struct quat q);

#ifdef __cplusplus
}
#endif

#endif	/* _LIBQUAT_H_ */
