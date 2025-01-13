/* this is the macro file for the adh model */

/* macros for adh */

#define ELEM2D_DOF3_COPY(f1, f2) \
        f1[0].x_eq = f2[0].x_eq;  \
        f1[1].x_eq = f2[1].x_eq;  \
        f1[2].x_eq = f2[2].x_eq;  \
        f1[0].y_eq = f2[0].y_eq;  \
        f1[1].y_eq = f2[1].y_eq;  \
        f1[2].y_eq = f2[2].y_eq;  \
        f1[0].c_eq = f2[0].c_eq;  \
        f1[1].c_eq = f2[1].c_eq;  \
        f1[2].c_eq = f2[2].c_eq;  \


#define LOCAL_INIT_VECT2D(local) local.x = 0.; local.y = 0.; 

#define VECT2D_MAG(vec) sqrt(vec.x * vec.x + vec.y * vec.y);
#define VECT2D_MAG_SAFE(vec) sqrt(vec.x * vec.x + vec.y * vec.y + NOT_QUITE_SMALL);
#define VECT3D_MAG(vec) sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
#define VECT3D_MAG_SAFE(vec) sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z + NOT_QUITE_SMALL);

#define ELEM1D_LOCAL_VECT2D_AVG(vect, ovg) \
        avg.x = (1./2.) * (vect[0].x + vect[1].x); \
        avg.y = (1./2.) * (vect[0].y + vect[1].y);

#define ELEM2D_LOCAL_VECT2D_AVG(vect, avg) \
        avg.x = (1./3.) * (vect[0].x + vect[1].x + vect[2].x); \
        avg.y = (1./3.) * (vect[0].y + vect[1].y + vect[2].y);

#define ELEM3D_LOCAL_VECT2D_AVG(vect, avg) \
        avg.x = (1./4.) * (vect[0].x + vect[1].x + vect[2].x + vect[3].x); \
        avg.y = (1./4.) * (vect[0].y + vect[1].y + vect[2].y + vect[3].y); 


#define ELEM3D_LOCAL_VECT_AVG(vect, avg) \
        avg.x = (1./4.) * (vect[0].x + vect[1].x + vect[2].x + vect[3].x); \
        avg.y = (1./4.) * (vect[0].y + vect[1].y + vect[2].y + vect[3].y); \
        avg.z = (1./4.) * (vect[0].z + vect[1].z + vect[2].z + vect[3].z);

#define ELEM3D_LOCAL_AVG(func, avg) \
        avg = (1./4.) * (func[0] + func[1] + func[2] + func[3]);

#define ELEM3D_LOCAL_VECT_SUM(vect, sum) \
        sum.x = vect[0].x + vect[1].x + vect[2].x + vect[3].x; \
        sum.y = vect[0].y + vect[1].y + vect[2].y + vect[3].y; \
        sum.z = vect[0].z + vect[1].z + vect[2].z + vect[3].z;

#define ELEM3D_LOCAL_SUM(func, sum) \
        sum = func[0] + func[1] + func[2] + func[3];

#define GRAD2D_FUNC(grad_shp, var, grad) \
        grad.x = var[0] * grad_shp[0].x + var[1] * grad_shp[1].x  + var[2] * grad_shp[2].x;\
        grad.y = var[0] * grad_shp[0].y + var[1] * grad_shp[1].y  + var[2] * grad_shp[2].y;

#define GRAD2D_VECT2D_X(grad_shp, var, grad) \
        grad.x = var[0].x * grad_shp[0].x + var[1].x * grad_shp[1].x  + var[2].x * grad_shp[2].x;\
        grad.y = var[0].x * grad_shp[0].y + var[1].x * grad_shp[1].y  + var[2].x * grad_shp[2].y;

#define GRAD2D_VECT2D_Y(grad_shp, var, grad) \
        grad.x = var[0].y * grad_shp[0].x + var[1].y * grad_shp[1].x  + var[2].y * grad_shp[2].x;\
        grad.y = var[0].y * grad_shp[0].y + var[1].y * grad_shp[1].y  + var[2].y * grad_shp[2].y;

#define VECT3D_MAG(vec) sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);

#define GRAD3D_FUNC(grad_shp, var, grad) \
        grad.x = var[0] * grad_shp[0].x + var[1] * grad_shp[1].x  + var[2] * grad_shp[2].x + var[3] * grad_shp[3].x; \
        grad.y = var[0] * grad_shp[0].y + var[1] * grad_shp[1].y  + var[2] * grad_shp[2].y + var[3] * grad_shp[3].y; \
        grad.z = var[0] * grad_shp[0].z + var[1] * grad_shp[1].z  + var[2] * grad_shp[2].z + var[3] * grad_shp[3].z; 

#define GRAD3D_VECT3D_X(grad_shp, var, grad) \
        grad.x = var[0].x * grad_shp[0].x + var[1].x * grad_shp[1].x  + var[2].x * grad_shp[2].x + var[3].x * grad_shp[3].x;\
        grad.y = var[0].x * grad_shp[0].y + var[1].x * grad_shp[1].y  + var[2].x * grad_shp[2].y + var[3].x * grad_shp[3].y;\
        grad.z = var[0].x * grad_shp[0].z + var[1].x * grad_shp[1].z  + var[2].x * grad_shp[2].z + var[3].x * grad_shp[3].z;

#define GRAD3D_VECT3D_Y(grad_shp, var, grad) \
        grad.x = var[0].y * grad_shp[0].x + var[1].y * grad_shp[1].x  + var[2].y * grad_shp[2].x + var[3].y * grad_shp[3].x;\
        grad.y = var[0].y * grad_shp[0].y + var[1].y * grad_shp[1].y  + var[2].y * grad_shp[2].y + var[3].y * grad_shp[3].y;\
        grad.z = var[0].y * grad_shp[0].z + var[1].y * grad_shp[1].z  + var[2].y * grad_shp[2].z + var[3].y * grad_shp[3].z;

#define GRAD3D_VECT3D_Z(grad_shp, var, grad) \
        grad.x = var[0].z * grad_shp[0].x + var[1].z * grad_shp[1].x  + var[2].z * grad_shp[2].x + var[3].z * grad_shp[3].x;\
        grad.y = var[0].z * grad_shp[0].y + var[1].z * grad_shp[1].y  + var[2].z * grad_shp[2].y + var[3].z * grad_shp[3].y;\
        grad.z = var[0].z * grad_shp[0].z + var[1].z * grad_shp[1].z  + var[2].z * grad_shp[2].z + var[3].z * grad_shp[3].z;

#define NUM_GET_TPOSITION(local1,local2,local3,scale,time1,time2) \
        local1=scale*(1.5*local2-0.5*local3)-(1.-scale)*time1*local2/time2; 

#define NUM_GET_TPOSITION_VECT3D(local1,local2,local3,scale,time1,time2) \
        local1.x = scale*(1.5*local2.x - 0.5*local3.x) - (1.-scale)*time1 * local2.x/time2; \
        local1.y = scale*(1.5*local2.y - 0.5*local3.y) - (1.-scale)*time1 * local2.y/time2; \
        local1.z = scale*(1.5*local2.z - 0.5*local3.z) - (1.-scale)*time1 * local2.z/time2;

#define NUM_GET_TPOSITION_VECT2D(local1,local2,local3,scale,time1,time2) \
        local1.x = scale*(1.5*local2.x - 0.5*local3.x) - (1.-scale)*time1 * local2.x/time2; \
        local1.y = scale*(1.5*local2.y - 0.5*local3.y) - (1.-scale)*time1 * local2.y/time2; 

#define NUM_GET_SECONDORDER(local1, local2, local3) \
        local1=(1.5*local2-0.5*local3);

#define ELEM3D_GET_VECT_TPOSITION(local1,local2,local3,scale,time1,time2) \
        local1[0].x=scale*(1.5*local2[0].x-0.5*local3[0].x)+(1.-scale)*time1*local2[0].x/time2; \
        local1[0].y=scale*(1.5*local2[0].y-0.5*local3[0].y)+(1.-scale)*time1*local2[0].y/time2; \
        local1[0].z=scale*(1.5*local2[0].z-0.5*local3[0].z)+(1.-scale)*time1*local2[0].z/time2; \
        local1[1].x=scale*(1.5*local2[1].x-0.5*local3[1].x)+(1.-scale)*time1*local2[1].x/time2; \
        local1[1].y=scale*(1.5*local2[1].y-0.5*local3[1].y)+(1.-scale)*time1*local2[1].y/time2; \
        local1[1].z=scale*(1.5*local2[1].z-0.5*local3[1].z)+(1.-scale)*time1*local2[1].z/time2; \
        local1[2].x=scale*(1.5*local2[2].x-0.5*local3[2].x)+(1.-scale)*time1*local2[2].x/time2; \
        local1[2].y=scale*(1.5*local2[2].y-0.5*local3[2].y)+(1.-scale)*time1*local2[2].y/time2; \
        local1[2].z=scale*(1.5*local2[2].z-0.5*local3[2].z)+(1.-scale)*time1*local2[2].z/time2; \
        local1[3].x=scale*(1.5*local2[3].x-0.5*local3[3].x)+(1.-scale)*time1*local2[3].x/time2; \
        local1[3].y=scale*(1.5*local2[3].y-0.5*local3[3].y)+(1.-scale)*time1*local2[3].y/time2; \
        local1[3].z=scale*(1.5*local2[3].z-0.5*local3[3].z)+(1.-scale)*time1*local2[3].z/time2;

#define ELEM3D_GET_VECT_SECONDORDER(local1,local2,local3) \
        local1[0].x=(1.5*local2[0].x-0.5*local3[0].x); \
        local1[0].y=(1.5*local2[0].y-0.5*local3[0].y); \
        local1[0].z=(1.5*local2[0].z-0.5*local3[0].z); \
        local1[1].x=(1.5*local2[1].x-0.5*local3[1].x); \
        local1[1].y=(1.5*local2[1].y-0.5*local3[1].y); \
        local1[1].z=(1.5*local2[1].z-0.5*local3[1].z); \
        local1[2].x=(1.5*local2[2].x-0.5*local3[2].x); \
        local1[2].y=(1.5*local2[2].y-0.5*local3[2].y); \
        local1[2].z=(1.5*local2[2].z-0.5*local3[2].z); \
        local1[3].x=(1.5*local2[3].x-0.5*local3[3].x); \
        local1[3].y=(1.5*local2[3].y-0.5*local3[3].y); \
        local1[3].z=(1.5*local2[3].z-0.5*local3[3].z);

#define ELEM3D_GET_TPOSITION(local1,local2,local3,scale,time1,time2) \
        local1[0]=scale*(1.5*local2[0]-0.5*local3[0])+(1.-scale)*time1*local2[0]/time2; \
        local1[1]=scale*(1.5*local2[1]-0.5*local3[1])+(1.-scale)*time1*local2[1]/time2; \
        local1[2]=scale*(1.5*local2[2]-0.5*local3[2])+(1.-scale)*time1*local2[2]/time2; \
        local1[3]=scale*(1.5*local2[3]-0.5*local3[3])+(1.-scale)*time1*local2[3]/time2;

#define ELEM3D_GET_TPOSITION_VECT(local1,local2,local3,scale,time1,time2) \
        local1[0].x=scale*(1.5*local2[0].x-0.5*local3[0].x)+(1.-scale)*time1*local2[0].x/time2; \
        local1[1].x=scale*(1.5*local2[1].x-0.5*local3[1].x)+(1.-scale)*time1*local2[1].x/time2; \
        local1[2].x=scale*(1.5*local2[2].x-0.5*local3[2].x)+(1.-scale)*time1*local2[2].x/time2; \
        local1[3].x=scale*(1.5*local2[3].x-0.5*local3[3].x)+(1.-scale)*time1*local2[3].x/time2; \
        local1[0].y=scale*(1.5*local2[0].y-0.5*local3[0].y)+(1.-scale)*time1*local2[0].y/time2; \
        local1[1].y=scale*(1.5*local2[1].y-0.5*local3[1].y)+(1.-scale)*time1*local2[1].y/time2; \
        local1[2].y=scale*(1.5*local2[2].y-0.5*local3[2].y)+(1.-scale)*time1*local2[2].y/time2; \
        local1[3].y=scale*(1.5*local2[3].y-0.5*local3[3].y)+(1.-scale)*time1*local2[3].y/time2; \
        local1[0].z=scale*(1.5*local2[0].z-0.5*local3[0].z)+(1.-scale)*time1*local2[0].z/time2; \
        local1[1].z=scale*(1.5*local2[1].z-0.5*local3[1].z)+(1.-scale)*time1*local2[1].z/time2; \
        local1[2].z=scale*(1.5*local2[2].z-0.5*local3[2].z)+(1.-scale)*time1*local2[2].z/time2; \
        local1[3].z=scale*(1.5*local2[3].z-0.5*local3[3].z)+(1.-scale)*time1*local2[3].z/time2; 

#define ELEM3D_GET_TPOSITION_VECT2D(local1,local2,local3,scale,time1,time2) \
        local1[0].x=scale*(1.5*local2[0].x-0.5*local3[0].x)+(1.-scale)*time1*local2[0].x/time2; \
        local1[1].x=scale*(1.5*local2[1].x-0.5*local3[1].x)+(1.-scale)*time1*local2[1].x/time2; \
        local1[2].x=scale*(1.5*local2[2].x-0.5*local3[2].x)+(1.-scale)*time1*local2[2].x/time2; \
        local1[3].x=scale*(1.5*local2[3].x-0.5*local3[3].x)+(1.-scale)*time1*local2[3].x/time2; \
        local1[0].y=scale*(1.5*local2[0].y-0.5*local3[0].y)+(1.-scale)*time1*local2[0].y/time2; \
        local1[1].y=scale*(1.5*local2[1].y-0.5*local3[1].y)+(1.-scale)*time1*local2[1].y/time2; \
        local1[2].y=scale*(1.5*local2[2].y-0.5*local3[2].y)+(1.-scale)*time1*local2[2].y/time2; \
        local1[3].y=scale*(1.5*local2[3].y-0.5*local3[3].y)+(1.-scale)*time1*local2[3].y/time2; 

#define ELEM3D_GET_SECONDORDER(local1,local2,local3) \
        local1[0]=(1.5*local2[0]-0.5*local3[0]); \
        local1[1]=(1.5*local2[1]-0.5*local3[1]); \
        local1[2]=(1.5*local2[2]-0.5*local3[2]); \
        local1[3]=(1.5*local2[3]-0.5*local3[3]);

#define ELEM2D_GET_SECONDORDER(local1,local2,local3) \
        local1[0]=(1.5*local2[0]-0.5*local3[0]); \
        local1[1]=(1.5*local2[1]-0.5*local3[1]); \
        local1[2]=(1.5*local2[2]-0.5*local3[2]);

#define ELEM2D_GET_TPOSITION(local1,local2,local3,scale,time1,time2) \
        local1[0]=scale*(1.5*local2[0]-0.5*local3[0])+(1.-scale)*time1*local2[0]/time2; \
        local1[1]=scale*(1.5*local2[1]-0.5*local3[1])+(1.-scale)*time1*local2[1]/time2; \
        local1[2]=scale*(1.5*local2[2]-0.5*local3[2])+(1.-scale)*time1*local2[2]/time2;

#define ELEM3D_GET_LOCAL(global, local, nodes) \
        local[0] = global[nodes[0]]; \
        local[1] = global[nodes[1]]; \
        local[2] = global[nodes[2]]; \
        local[3] = global[nodes[3]];

#define ELEM2D_GET_LOCAL(global, local, nodes) \
        local[0] = global[nodes[0]]; \
        local[1] = global[nodes[1]]; \
        local[2] = global[nodes[2]];

#define ELEM1D_GET_LOCAL(global, local, nodes) \
        local[0] = global[nodes[0]]; \
        local[1] = global[nodes[1]];

#define ELEM3D_GET_LOCAL_VECT(global, local, nodes) \
        local[0].x = global[nodes[0]].x; \
        local[0].y = global[nodes[0]].y; \
        local[0].z = global[nodes[0]].z; \
        local[1].x = global[nodes[1]].x; \
        local[1].y = global[nodes[1]].y; \
        local[1].z = global[nodes[1]].z; \
        local[2].x = global[nodes[2]].x; \
        local[2].y = global[nodes[2]].y; \
        local[2].z = global[nodes[2]].z; \
        local[3].x = global[nodes[3]].x; \
        local[3].y = global[nodes[3]].y; \
        local[3].z = global[nodes[3]].z;

#define ELEM3D_GET_LOCAL_VECT2D(global, local, nodes) \
        local[0].x = global[nodes[0]].x; \
        local[0].y = global[nodes[0]].y; \
        local[1].x = global[nodes[1]].x; \
        local[1].y = global[nodes[1]].y; \
        local[2].x = global[nodes[2]].x; \
        local[2].y = global[nodes[2]].y; \
        local[3].x = global[nodes[3]].x; \
        local[3].y = global[nodes[3]].y;

#define ELEM2D_GET_LOCAL_VECT(global, local, nodes) \
        local[0].x = global[nodes[0]].x; \
        local[0].y = global[nodes[0]].y; \
        local[0].z = global[nodes[0]].z; \
        local[1].x = global[nodes[1]].x; \
        local[1].y = global[nodes[1]].y; \
        local[1].z = global[nodes[1]].z; \
        local[2].x = global[nodes[2]].x; \
        local[2].y = global[nodes[2]].y; \
        local[2].z = global[nodes[2]].z;

#define ELEM1D_GET_LOCAL_VECT(global, local, nodes) \
        local[0].x = global[nodes[0]].x; \
        local[0].y = global[nodes[0]].y; \
        local[0].z = global[nodes[0]].z; \
        local[1].x = global[nodes[1]].x; \
        local[1].y = global[nodes[1]].y; \
        local[1].z = global[nodes[1]].z;

#define ELEM2D_GET_LOCAL_TENSOR2D(global, local, nodes) \
        local[0].xx = global[nodes[0]].xx; \
        local[0].yy = global[nodes[0]].yy; \
        local[0].xy = global[nodes[0]].xy; \
        local[1].xx = global[nodes[1]].xx; \
        local[1].yy = global[nodes[1]].yy; \
        local[1].xy = global[nodes[1]].xy; \
        local[2].xx = global[nodes[2]].xx; \
        local[2].yy = global[nodes[2]].yy; \
        local[2].xy = global[nodes[2]].xy;

#define ELEM2D_GET_LOCAL_VECT2D(global, local, nodes) \
        local[0].x = global[nodes[0]].x; \
        local[0].y = global[nodes[0]].y; \
        local[1].x = global[nodes[1]].x; \
        local[1].y = global[nodes[1]].y; \
        local[2].x = global[nodes[2]].x; \
        local[2].y = global[nodes[2]].y;

#define ELEM1D_GET_LOCAL_VECT2D(global, local, nodes) \
        local[0].x = global[nodes[0]].x; \
        local[0].y = global[nodes[0]].y; \
        local[1].x = global[nodes[1]].x; \
        local[1].y = global[nodes[1]].y;

#define ELEM3D_GET_ELEVATION(local_z, nodes) \
        local_z[0] = node[nodes[0]].z; \
        local_z[1] = node[nodes[1]].z; \
        local_z[2] = node[nodes[2]].z; \
        local_z[3] = node[nodes[3]].z;

#define ELEM2D_GET_ELEVATION(local_z, nodes) \
        local_z[0] = node[nodes[0]].z; \
        local_z[1] = node[nodes[1]].z; \
        local_z[2] = node[nodes[2]].z;

#define ELEM1D_GET_ELEVATION(local_z, nodes) \
        local_z[0] = node[nodes[0]].z; \
        local_z[1] = node[nodes[1]].z;

#define ELEM3D_GET_SUM(local) (local[0]+local[1]+local[2]+local[3])

#define ELEM3D_GET_AVG(global, nodes) (0.25*(global[nodes[0]]+global[nodes[1]]+global[nodes[2]]+global[nodes[3]]))

#define ELEM2D_GET_AVG(global, nodes) ((global[nodes[0]]+global[nodes[1]]+global[nodes[2]])/3.0)

#define ELEM1D_GET_AVG(global, nodes) (0.5*(global[nodes[0]]+global[nodes[1]]))

#define ELEM3D_ASSMB_M_1DOF(global, local, nodes) \
        global[nodes[0]] -= local[0]; \
        global[nodes[1]] -= local[1]; \
        global[nodes[2]] -= local[2]; \
        global[nodes[3]] -= local[3];

#define ELEM2D_ASSMB_M_1DOF(global, local, nodes) \
        global[nodes[0]] -= local[0]; \
        global[nodes[1]] -= local[1]; \
        global[nodes[2]] -= local[2];

#define ELEM1D_ASSMB_M_1DOF(global, local, nodes) \
        global[nodes[0]] -= local[0]; \
        global[nodes[1]] -= local[1];

#define ELEM3D_LOCAL_COPY(local1, local2) \
        local2[0] = local1[0]; \
        local2[1] = local1[1]; \
        local2[2] = local1[2]; \
        local2[3] = local1[3];

#define ELEM2D_LOCAL_COPY(local1, local2) \
        local2[0] = local1[0]; \
        local2[1] = local1[1]; \
        local2[2] = local1[2];

#define ELEM1D_LOCAL_COPY(local1, local2) \
        local2[0] = local1[0]; \
        local2[1] = local1[1];

#define ELEM3D_GET_VECT_SUM(local1, local2) \
        local2.x = local1[0].x+local1[1].x+local1[2].x+local1[3].x; \
        local2.y = local1[0].y+local1[1].y+local1[2].y+local1[3].y; \
        local2.z = local1[0].z+local1[1].z+local1[2].z+local1[3].z;

#define ELEM2D_GET_VECT_SUM(local1, local2) \
        local2.x = local1[0].x+local1[1].x+local1[2].x; \
        local2.y = local1[0].y+local1[1].y+local1[2].y; \
        local2.z = local1[0].z+local1[1].z+local1[2].z;

#define ELEM2D_GET_VECT2D_SUM(local1, local2) \
        local2.x = local1[0].x+local1[1].x+local1[2].x; \
        local2.y = local1[0].y+local1[1].y+local1[2].y;

#define ELEM3D_VECT_DIFF(result,local1, local2) \
        result[0].x = local1[0].x-local2[0].x; \
        result[1].x = local1[1].x-local2[1].x; \
        result[2].x = local1[2].x-local2[2].x; \
        result[3].x = local1[3].x-local2[3].x; \
        result[0].y = local1[0].y-local2[0].y; \
        result[1].y = local1[1].y-local2[1].y; \
        result[2].y = local1[2].y-local2[2].y; \
        result[3].y = local1[3].y-local2[3].y; \
        result[0].z = local1[0].z-local2[0].z; \
        result[1].z = local1[1].z-local2[1].z; \
        result[2].z = local1[2].z-local2[2].z; \
        result[3].z = local1[3].z-local2[3].z;

#define ELEM2D_VECT_DIFF(result,local1, local2) \
        result[0].x = local1[0].x-local2[0].x; \
        result[1].x = local1[1].x-local2[1].x; \
        result[2].x = local1[2].x-local2[2].x; \
        result[0].y = local1[0].y-local2[0].y; \
        result[1].y = local1[1].y-local2[1].y; \
        result[2].y = local1[2].y-local2[2].y; \
        result[0].z = local1[0].z-local2[0].z; \
        result[1].z = local1[1].z-local2[1].z; \
        result[2].z = local1[2].z-local2[2].z;

#define ELEM2D_VECT2D_SCALE(local,scale) \
        local[0].x *= scale; \
        local[0].y *= scale; \
        local[1].x *= scale; \
        local[1].y *= scale; \
        local[2].x *= scale; \
        local[2].y *= scale;

#define ELEM3D_VECT_SCALE(local,scale) \
        local[0].x *= scale; \
        local[0].y *= scale; \
        local[0].z *= scale; \
        local[1].x *= scale; \
        local[1].y *= scale; \
        local[1].z *= scale; \
        local[2].x *= scale; \
        local[2].y *= scale; \
        local[2].z *= scale; \
        local[3].x *= scale; \
        local[3].y *= scale; \
        local[3].z *= scale;

#define ELEM2D_VECT_SCALE(local,scale) \
        local[0].x *= scale; \
        local[0].y *= scale; \
        local[0].z *= scale; \
        local[1].x *= scale; \
        local[1].y *= scale; \
        local[1].z *= scale; \
        local[2].x *= scale; \
        local[2].y *= scale; \
        local[2].z *= scale;

#define ELEM1D_VECT_SCALE(local,scale) \
        local[0].x *= scale; \
        local[0].y *= scale; \
        local[0].z *= scale; \
        local[1].x *= scale; \
        local[1].y *= scale; \
        local[1].z *= scale;

#define ELEM3D_SCALE(local,scale) \
        local[0] *= scale; \
        local[1] *= scale; \
        local[2] *= scale; \
        local[3] *= scale;

#define ELEM2D_SCALE(local,scale) \
        local[0] *= scale; \
        local[1] *= scale; \
        local[2] *= scale;

#define ELEM1D_SCALE(local,scale) \
        local[0] *= scale; \
        local[1] *= scale;

#define ELEM3D_LOCAL_VECT_COPY(local1, local2) \
        local2[0].x = local1[0].x; \
        local2[1].x = local1[1].x; \
        local2[2].x = local1[2].x; \
        local2[3].x = local1[3].x; \
        local2[0].y = local1[0].y; \
        local2[1].y = local1[1].y; \
        local2[2].y = local1[2].y; \
        local2[3].y = local1[3].y; \
        local2[0].z = local1[0].z; \
        local2[1].z = local1[1].z; \
        local2[2].z = local1[2].z; \
        local2[3].z = local1[3].z;

#define ELEM2D_LOCAL_VECT2D_COPY(local1, local2) \
        local2[0].x = local1[0].x; \
        local2[1].x = local1[1].x; \
        local2[2].x = local1[2].x; \
        local2[0].y = local1[0].y; \
        local2[1].y = local1[1].y; \
        local2[2].y = local1[2].y;

#define ELEM1D_LOCAL_VECT2D_COPY(local1, local2) \
        local2[0].x = local1[0].x; \
        local2[1].x = local1[1].x; \
        local2[0].y = local1[0].y; \
        local2[1].y = local1[1].y;

#define ELEM2D_LOCAL_VECT_COPY(local1, local2) \
        local2[0].x = local1[0].x; \
        local2[1].x = local1[1].x; \
        local2[2].x = local1[2].x; \
        local2[0].y = local1[0].y; \
        local2[1].y = local1[1].y; \
        local2[2].y = local1[2].y; \
        local2[0].z = local1[0].z; \
        local2[1].z = local1[1].z; \
        local2[2].z = local1[2].z;

#define ELEM3D_LOCAL_DIFF(result, local1, local2) \
        result[0] = local1[0]-local2[0]; \
        result[1] = local1[1]-local2[1]; \
        result[2] = local1[2]-local2[2]; \
        result[3] = local1[3]-local2[3];

#define ELEM2D_LOCAL_DIFF(result, local1, local2) \
        result[0] = local1[0]-local2[0]; \
        result[1] = local1[1]-local2[1]; \
        result[2] = local1[2]-local2[2];

#define ELEM1D_LOCAL_DIFF(result, local1, local2) \
        result[0] = local1[0]-local2[0]; \
        result[1] = local1[1]-local2[1];

#define ELEM2D_LOCAL_ADD(result, local1, local2)\
        result[0] = local1[0]+local2[0]; \
        result[1] = local1[1]+local2[1]; \
        result[2] = local1[2]+local2[2];

#define ELEM1D_LOCAL_ADD(result, local1, local2) \
        result[0] = local1[0]+local2[0]; \
        result[1] = local1[1]+local2[1];


#define ELEM3D_INIT_MAT_1DOF(emat) \
        emat[ 0] = 1.0; emat[ 1] = 0.0; emat[ 2] = 0.0; emat[ 3] = 0.0; \
        emat[ 4] = 0.0; emat[ 5] = 1.0; emat[ 6] = 0.0; emat[ 7] = 0.0; \
        emat[ 8] = 0.0; emat[ 9] = 0.0; emat[10] = 1.0; emat[11] = 0.0; \
        emat[12] = 0.0; emat[13] = 0.0; emat[14] = 0.0; emat[15] = 1.0;

#define ELEM3D_ZERO_MAT_1DOF(emat) \
        emat[ 0] = 0.0; emat[ 1] = 0.0; emat[ 2] = 0.0; emat[ 3] = 0.0; \
        emat[ 4] = 0.0; emat[ 5] = 0.0; emat[ 6] = 0.0; emat[ 7] = 0.0; \
        emat[ 8] = 0.0; emat[ 9] = 0.0; emat[10] = 0.0; emat[11] = 0.0; \
        emat[12] = 0.0; emat[13] = 0.0; emat[14] = 0.0; emat[15] = 0.0;

#define ELEM2D_INIT_MAT_1DOF(emat) \
        emat[ 0] = 1.0; emat[ 1] = 0.0; emat[ 2] = 0.0; \
        emat[ 3] = 0.0; emat[ 4] = 1.0; emat[ 5] = 0.0; \
        emat[ 6] = 0.0; emat[ 7] = 0.0; emat[ 8] = 1.0;

#define ELEM2D_ZERO_MAT_1DOF(emat) \
        emat[ 0] = 0.0; emat[ 1] = 0.0; emat[ 2] = 0.0; \
        emat[ 3] = 0.0; emat[ 4] = 0.0; emat[ 5] = 0.0; \
        emat[ 6] = 0.0; emat[ 7] = 0.0; emat[ 8] = 0.0;

#define ELEM2D_INIT_MAT_3DOF(emat) \
        emat[ 0].x_eq = 1.0; emat[ 1].x_eq = 0.0; emat[ 2].x_eq = 0.0; \
        emat[ 3].x_eq = 0.0; emat[ 4].x_eq = 1.0; emat[ 5].x_eq = 0.0; \
        emat[ 6].x_eq = 0.0; emat[ 7].x_eq = 0.0; emat[ 8].x_eq = 1.0; \
        emat[ 0].y_eq = 1.0; emat[ 1].y_eq = 0.0; emat[ 2].y_eq = 0.0; \
        emat[ 3].y_eq = 0.0; emat[ 4].y_eq = 1.0; emat[ 5].y_eq = 0.0; \
        emat[ 6].y_eq = 0.0; emat[ 7].y_eq = 0.0; emat[ 8].y_eq = 1.0; \
        emat[ 0].c_eq = 1.0; emat[ 1].c_eq = 0.0; emat[ 2].c_eq = 0.0; \
        emat[ 3].c_eq = 0.0; emat[ 4].c_eq = 1.0; emat[ 5].c_eq = 0.0; \
        emat[ 6].c_eq = 0.0; emat[ 7].c_eq = 0.0; emat[ 8].c_eq = 1.0;

#define ELEM3D_ZERO_MAT_3DOF(emat) \
        emat[0].x_eq = 0.; emat[1].x_eq=0.; emat[2].x_eq=0.; emat[3].x_eq=0.; \
        emat[4].x_eq = 0.; emat[5].x_eq=0.; emat[6].x_eq=0.; emat[7].x_eq=0.; \
        emat[8].x_eq = 0.; emat[9].x_eq=0.; emat[10].x_eq=0.; emat[11].x_eq=0.; \
        emat[12].x_eq = 0.; emat[13].x_eq=0.; emat[14].x_eq=0.; emat[15].x_eq=0.;\
        emat[0].y_eq = 0.; emat[1].y_eq=0.; emat[2].y_eq=0.; emat[3].y_eq=0.; \
        emat[4].y_eq = 0.; emat[5].y_eq=0.; emat[6].y_eq=0.; emat[7].y_eq=0.; \
        emat[8].y_eq = 0.; emat[9].y_eq=0.; emat[10].y_eq=0.; emat[11].y_eq=0.; \
        emat[12].y_eq = 0.; emat[13].y_eq=0.; emat[14].y_eq=0.; emat[15].y_eq=0.;\
        emat[0].c_eq = 0.; emat[1].c_eq=0.; emat[2].c_eq=0.; emat[3].c_eq=0.; \
        emat[4].c_eq = 0.; emat[5].c_eq=0.; emat[6].c_eq=0.; emat[7].c_eq=0.; \
        emat[8].c_eq = 0.; emat[9].c_eq=0.; emat[10].c_eq=0.; emat[11].c_eq=0.; \
        emat[12].c_eq = 0.; emat[13].c_eq=0.; emat[14].c_eq=0.; emat[15].c_eq=0.;

#define ELEM2D_ZERO_MAT_3DOF(emat) \
        emat[0].x_eq = 0.; emat[1].x_eq = 0.; emat[2].x_eq = 0.;  \
        emat[3].x_eq = 0.; emat[4].x_eq = 0.; emat[5].x_eq = 0.;  \
        emat[6].x_eq = 0.; emat[7].x_eq = 0.; emat[8].x_eq = 0.;  \
        emat[0].y_eq = 0.; emat[1].y_eq = 0.; emat[2].y_eq = 0.;  \
        emat[3].y_eq = 0.; emat[4].y_eq = 0.; emat[5].y_eq = 0.;  \
        emat[6].y_eq = 0.; emat[7].y_eq = 0.; emat[8].y_eq = 0.;  \
        emat[0].c_eq = 0.; emat[1].c_eq = 0.; emat[2].c_eq = 0.;  \
        emat[3].c_eq = 0.; emat[4].c_eq = 0.; emat[5].c_eq = 0.;  \
        emat[6].c_eq = 0.; emat[7].c_eq = 0.; emat[8].c_eq = 0.;  

#define ELEM1D_INIT_MAT_3DOF(emat) \
        emat[ 0].x_eq = 1.0; emat[ 1].x_eq = 0.0; \
        emat[ 2].x_eq = 0.0; emat[ 3].x_eq = 1.0; \
        emat[ 0].y_eq = 1.0; emat[ 1].y_eq = 0.0; \
        emat[ 2].y_eq = 0.0; emat[ 3].y_eq = 1.0; \
        emat[ 0].c_eq = 1.0; emat[ 1].c_eq = 0.0; \
        emat[ 2].c_eq = 0.0; emat[ 3].c_eq = 1.0;

#define ELEM1D_INIT_MAT_1DOF(emat) \
        emat[ 0] = 1.0; emat[ 1] = 0.0; \
        emat[ 2] = 0.0; emat[ 3] = 1.0;

#define ELEM3D_DERIV_1DOF(emat, local1, local2, index, diff_ep) \
        emat[   index] = (local1[0]-local2[0])/diff_ep; \
        emat[ 4+index] = (local1[1]-local2[1])/diff_ep; \
        emat[ 8+index] = (local1[2]-local2[2])/diff_ep; \
        emat[12+index] = (local1[3]-local2[3])/diff_ep;

#define ELEM3D_DERIV_3DOF(emat, local1, local2, index, diff_ep) \
        emat[   index].x_eq = (local1[0].x_eq-local2[0].x_eq)/diff_ep; \
        emat[ 4+index].x_eq = (local1[1].x_eq-local2[1].x_eq)/diff_ep; \
        emat[ 8+index].x_eq = (local1[2].x_eq-local2[2].x_eq)/diff_ep; \
        emat[12+index].x_eq = (local1[3].x_eq-local2[3].x_eq)/diff_ep; \
        emat[   index].y_eq = (local1[0].y_eq-local2[0].y_eq)/diff_ep; \
        emat[ 4+index].y_eq = (local1[1].y_eq-local2[1].y_eq)/diff_ep; \
        emat[ 8+index].y_eq = (local1[2].y_eq-local2[2].y_eq)/diff_ep; \
        emat[12+index].y_eq = (local1[3].y_eq-local2[3].y_eq)/diff_ep; \
        emat[   index].c_eq = (local1[0].c_eq-local2[0].c_eq)/diff_ep; \
        emat[ 4+index].c_eq = (local1[1].c_eq-local2[1].c_eq)/diff_ep; \
        emat[ 8+index].c_eq = (local1[2].c_eq-local2[2].c_eq)/diff_ep; \
        emat[12+index].c_eq = (local1[3].c_eq-local2[3].c_eq)/diff_ep; 

#define ELEM3D_DERIV_3DOF_2(emat, local1, local2, index, diff_ep_inv) \
        emat[   index].x_eq = diff_ep_inv * (local1[0].x_eq-local2[0].x_eq); \
        emat[ 4+index].x_eq = diff_ep_inv * (local1[1].x_eq-local2[1].x_eq); \
        emat[ 8+index].x_eq = diff_ep_inv * (local1[2].x_eq-local2[2].x_eq); \
        emat[12+index].x_eq = diff_ep_inv * (local1[3].x_eq-local2[3].x_eq); \
        emat[   index].y_eq = diff_ep_inv * (local1[0].y_eq-local2[0].y_eq); \
        emat[ 4+index].y_eq = diff_ep_inv * (local1[1].y_eq-local2[1].y_eq); \
        emat[ 8+index].y_eq = diff_ep_inv * (local1[2].y_eq-local2[2].y_eq); \
        emat[12+index].y_eq = diff_ep_inv * (local1[3].y_eq-local2[3].y_eq); \
        emat[   index].c_eq = diff_ep_inv * (local1[0].c_eq-local2[0].c_eq); \
        emat[ 4+index].c_eq = diff_ep_inv * (local1[1].c_eq-local2[1].c_eq); \
        emat[ 8+index].c_eq = diff_ep_inv * (local1[2].c_eq-local2[2].c_eq); \
        emat[12+index].c_eq = diff_ep_inv * (local1[3].c_eq-local2[3].c_eq);

#define ELEM2D_DERIV_1DOF(emat, local1, local2, index, diff_ep) \
        emat[   index] = (local1[0]-local2[0])/diff_ep; \
        emat[ 3+index] = (local1[1]-local2[1])/diff_ep; \
        emat[ 6+index] = (local1[2]-local2[2])/diff_ep;

#define ELEM2D_DERIV_3DOF(emat, local1, local2, index, diff_ep) \
        emat[   index].x_eq = (local1[0].x_eq-local2[0].x_eq)/diff_ep; \
        emat[ 3+index].x_eq = (local1[1].x_eq-local2[1].x_eq)/diff_ep; \
        emat[ 6+index].x_eq = (local1[2].x_eq-local2[2].x_eq)/diff_ep; \
        emat[   index].y_eq = (local1[0].y_eq-local2[0].y_eq)/diff_ep; \
        emat[ 3+index].y_eq = (local1[1].y_eq-local2[1].y_eq)/diff_ep; \
        emat[ 6+index].y_eq = (local1[2].y_eq-local2[2].y_eq)/diff_ep; \
        emat[   index].c_eq = (local1[0].c_eq-local2[0].c_eq)/diff_ep; \
        emat[ 3+index].c_eq = (local1[1].c_eq-local2[1].c_eq)/diff_ep; \
        emat[ 6+index].c_eq = (local1[2].c_eq-local2[2].c_eq)/diff_ep;

#define ELEM1D_DERIV_3DOF(emat, local1, local2, index, diff_ep) \
        emat[   index].x_eq = (local1[0].x_eq-local2[0].x_eq)/diff_ep; \
        emat[ 2+index].x_eq = (local1[1].x_eq-local2[1].x_eq)/diff_ep; \
        emat[   index].y_eq = (local1[0].y_eq-local2[0].y_eq)/diff_ep; \
        emat[ 2+index].y_eq = (local1[1].y_eq-local2[1].y_eq)/diff_ep; \
        emat[   index].c_eq = (local1[0].c_eq-local2[0].c_eq)/diff_ep; \
        emat[ 2+index].c_eq = (local1[1].c_eq-local2[1].c_eq)/diff_ep;

#define ELEM1D_DERIV_1DOF(emat, local1, local2, index, diff_ep) \
        emat[   index] = (local1[0]-local2[0])/diff_ep; \
        emat[ 2+index] = (local1[1]-local2[1])/diff_ep;

#define ELEM1D_LOCAL_INIT_DBL(local) local[0] = 0.0; local[1] = 0.0;

#define ELEM2D_LOCAL_INIT_DBL(local) local[0] = 0.0; local[1] = 0.0; local[2] = 0.0;

#define ELEM2D_LOCAL_INIT_INT(local) local[0] = 0; local[1] = 0; local[2] = 0;

#define ELEM3D_LOCAL_INIT_INT(local) local[0] = 0; local[1] = 0; local[2] = 0; local[3] = 0;

#define ELEM3D_LOCAL_INIT_DBL(local) local[0] = 0.0; local[1] = 0.0; local[2] = 0.0; local[3] = 0.0;

#define ELEM3D_LOCAL_QUAD_INIT_DBL(local) local[0] = 0.0; local[1] = 0.0; local[2] = 0.0; local[3] = 0.0; \
        local[4] = 0.0; local[5] = 0.0; local[6] = 0.0; local[7] = 0.0; local[8] = 0.0; local[9] = 0.0;

#define ELEM1D_LOCAL_INIT_IDNT_DBL(local) local[0]=1.; local[1]=1.;

#define ELEM2D_LOCAL_INIT_IDNT_DBL(local) local[0]=1.; local[1]=1.; local[2]=1.;

#define ELEM3D_LOCAL_INIT_IDNT_DBL(local) local[0]=1.; local[1]=1.; local[2]=1.; local[3]=1.;

#define ELEM3D_LOCAL_INIT_4DOF(local) local[0].x_eq=0.; local[0].y_eq=0.; local[0].z_eq=0.; local[0].c_eq=0.; \
        local[1].x_eq=0.; local[1].y_eq=0.; local[1].z_eq=0.; local[1].c_eq=0.; \
        local[2].x_eq=0.; local[2].y_eq=0.; local[2].z_eq=0.; local[2].c_eq=0.; \
        local[3].x_eq=0.; local[3].y_eq=0.; local[3].z_eq=0.; local[3].c_eq=0.;

#define ELEM2D_LOCAL_INIT_3DOF(local) local[0].x_eq=0.; local[0].y_eq=0.; local[0].c_eq=0.; \
        local[1].x_eq=0.; local[1].y_eq=0.; local[1].c_eq=0.; \
        local[2].x_eq=0.; local[2].y_eq=0.; local[2].c_eq=0.;

#define ELEM1D_LOCAL_INIT_3DOF(local) local[0].x_eq=0.; local[0].y_eq=0.; local[0].c_eq=0.; \
        local[1].x_eq=0.; local[1].y_eq=0.; local[1].c_eq=0.;

#define ELEM3D_LOCAL_INIT_2DOF(local) local[0].x_eq=0.; local[0].y_eq=0.; \
        local[1].x_eq=0.; local[1].y_eq=0.; \
        local[2].x_eq=0.; local[2].y_eq=0.; \
        local[3].x_eq=0.; local[3].y_eq=0.;

#define ELEM3D_LOCAL_INIT_3DOF(local) local[0].x_eq=0.; local[0].y_eq=0.; local[0].c_eq=0.; \
        local[1].x_eq=0.; local[1].y_eq=0.; local[1].c_eq=0.; \
        local[2].x_eq=0.; local[2].y_eq=0.; local[2].c_eq=0.; \
        local[3].x_eq=0.; local[3].y_eq=0.; local[3].c_eq=0.;

#define ELEM3D_LOCAL_INIT_VECT(local) local[0].x=0.; local[0].y=0.; local[0].z=0.; \
        local[1].x=0.; local[1].y=0.; local[1].z=0.; \
        local[2].x=0.; local[2].y=0.; local[2].z=0.; \
        local[3].x=0.; local[3].y=0.; local[3].z=0.;

#define ELEM3D_LOCAL_INIT_VECT2D(local) local[0].x=0.; local[0].y=0.; \
        local[1].x=0.; local[1].y=0.; \
        local[2].x=0.; local[2].y=0.; \
        local[3].x=0.; local[3].y=0.; 

#define ELEM2D_LOCAL_INIT_VECT(local) local[0].x=0.; local[0].y=0.; local[0].z=0.; \
        local[1].x=0.; local[1].y=0.; local[1].z=0.; \
        local[2].x=0.; local[2].y=0.; local[2].z=0.;

#define ELEM2D_LOCAL_INIT_VECT2D(local) local[0].x = 0.; local[0].y = 0.; \
        local[1].x = 0.; local[1].y = 0.; local[2].x = 0.; local[2].y = 0.;

#define ELEM2D_LOCAL_INIT_TENSOR2D(local) local[0].xx = 0.; local[0].xy = 0.; local[0].yy = 0.; \
        local[1].xx = 0.; local[1].xy = 0.; local[1].yy = 0.; \
        local[2].xx = 0.; local[2].xy = 0.; local[2].yy = 0.;

#define ELEM2D_INTGRT_DBLDBL(local1,local2) (local1[0]*local2[0]+local1[1]*local2[1]+local1[2]*local2[2]+ \
          (local1[0]+local1[1]+local1[2])*(local2[0]+local2[1]+local2[2]))/12.

/* in NUM_DIFF_EPSILON we scale by nodal variable */
#define NUM_DIFF_EPSILON(diff_ep, diff_ep2, var, perturbation) \
        diff_ep = fabs(var)*(1.+perturbation); \
        diff_ep -= fabs(var); \
        if(diff_ep <perturbation ) diff_ep=perturbation; \
        diff_ep2 = 2.0*diff_ep;

/* NUM_DIFF_EPSILON_GENERAL we scale by nodla variable and user specified fraction  */
#define NUM_DIFF_EPSILON_GENERAL(diff_ep, diff_ep2, var1, var2) \
        diff_ep = fabs(var1)*(1.+var2); \
        diff_ep -= fabs(var1); \
        if(diff_ep < var2 ) diff_ep=var2; \
        diff_ep2 = 2.0*diff_ep;

#define MAX_ON_ELEMENT(local,var); \
        var = fabs(local[0]); \
                if(fabs(local[1])>var)var=fabs(local[1]); \
                if(fabs(local[2])>var)var=fabs(local[2]); \
                if(fabs(local[3])>var)var=fabs(local[3]);

#define MAX_ON_FACE(local,var); \
        var = fabs(local[0]); \
                if(fabs(local[1])>var)var=fabs(local[1]); \
                if(fabs(local[2])>var)var=fabs(local[2]); 

#define VT_3D_TCOPY(tens, target) \
        target.xx = tens.xx; \
        target.yy = tens.yy; \
        target.zz = tens.zz; \
        target.xy = tens.xy; \
        target.xz = tens.xz; \
        target.yz = tens.yz;

#define VT_3D_VCOPY(vect, target) \
        target.x = vect.x; \
        target.y = vect.y; \
        target.z = vect.z;

#define VT_3D_INIT(vect) \
        vect.x = DZERO; \
        vect.y = DZERO; \
        vect.z = DZERO;

#define VT_1D_DOT(v1, v2) (v1*v2)

#define VT_2D_DOT(v1, v2) (v1.x*v2.x+v1.y*v2.y)

#define VT_3D_DOT(v1, v2) (v1.x*v2.x+v1.y*v2.y+v1.z*v2.z)

#define VT_3D_TENS_VECT_PROD(result, tens, vect) \
        result.x = tens.xx*vect.x+tens.xy*vect.y+tens.xz*vect.z; \
        result.y = tens.xy*vect.x+tens.yy*vect.y+tens.yz*vect.z; \
        result.z = tens.xz*vect.x+tens.yz*vect.y+tens.zz*vect.z;

#define VT_2D_TENS2D_VECT2D_PROD(result, tens, vect) \
        result.x = tens.xx*vect.x+tens.xy*vect.y; \
        result.y = tens.xy*vect.x+tens.yy*vect.y;

#define VT_2D_TENS2D_VECT2D_PROD_AI(result, tens, vect, vdir) \
        result.x = tens.xx*vdir.x*(vdir.x*vect.x+vdir.y*vect.y) + tens.xy*vect.x; \
        result.y = tens.yy*vdir.y*(vdir.x*vect.x+vdir.y*vect.y) + tens.xy*vect.y;

#define VT_3D_TSCALE(tens, scale) \
        tens.xx *= scale; \
        tens.xy *= scale; \
        tens.xz *= scale; \
        tens.yy *= scale; \
        tens.yz *= scale; \
        tens.zz *= scale;

#define VT_3D_VSCALE(vect, scale) \
        vect.x *= scale; \
        vect.y *= scale; \
        vect.z *= scale; \

#define VT_3D_TSIGN(tens) \
        tens.xx = -tens.xx; \
        tens.xy = -tens.xy; \
        tens.xz = -tens.xz; \
        tens.yy = -tens.yy; \
        tens.yz = -tens.yz; \
        tens.zz = -tens.zz;

#define VT_3D_VSIGN(vect) \
        vect.x = -vect.x; \
        vect.y = -vect.y; \
        vect.z = -vect.z; \

#define FIXED_POS(loc1,loc2,weight,new_loc)\
{\
        new_loc = loc1 + weight*(loc2-loc1);\
}

#define DIR_LENGTH(loc1,loc2) (loc1-loc2)

#define DIR_LENGTH_VEC(v0,v1,v2,vec99) \
{ \
        vec99[0].x=DIR_LENGTH(v1.x,v0.x); \
        vec99[0].y=DIR_LENGTH(v1.y,v0.y); \
        vec99[1].x=DIR_LENGTH(v2.x,v0.x); \
        vec99[1].y=DIR_LENGTH(v2.y,v0.y); \
}

#define TRI_AREA(v0,v1,v2,area) \
{ \
        SVECT2D vec[2]; \
        DIR_LENGTH_VEC(v0,v1,v2,vec); \
        area = (vec[0].x*vec[1].y-vec[0].y*vec[1].x)/2.0; \
}

#define GRAD_SHAPE(v0,v1,v2,area99,grad99) \
{ \
        grad99[0].x = 0.5*(v1.y-v2.y)/area99; \
        grad99[1].x = 0.5*(v2.y-v0.y)/area99; \
        grad99[2].x = 0.5*(v0.y-v1.y)/area99; \
        grad99[0].y = 0.5*(v2.x-v1.x)/area99; \
        grad99[1].y = 0.5*(v0.x-v2.x)/area99; \
        grad99[2].y = 0.5*(v1.x-v0.x)/area99; \
}

#define RE_DISTRIBUTE(x99,new_x99,new_elem_rhs99,only99) \
{ \
        int     i99, j99; \
        double xsi99, l99[2]; \
        SVECT2D vec[2]; \
        for(i99=0;i99<3;i99++) { \
                if(i99==only99) continue; \
                DIR_LENGTH_VEC(x99[only99],new_x99[i99],x99[i99],vec); \
                for(j99=0;j99<2;j99++) \
                        l99[j99] = vec[j99].x*vec[j99].x + vec[j99].y*vec[j99].y; \
                xsi99 = sqrt(l99[0]/l99[1]); \
                new_elem_rhs99[only99*3] += (1.0-xsi99) * new_elem_rhs99[i99*3]; \
                new_elem_rhs99[only99*3+1] += (1.0-xsi99) * new_elem_rhs99[i99*3+1]; \
                new_elem_rhs99[only99*3+2] += (1.0-xsi99) * new_elem_rhs99[i99*3+2]; \
                new_elem_rhs99[i99*3] *= xsi99; \
                new_elem_rhs99[i99*3+1] *= xsi99; \
                new_elem_rhs99[i99*3+2] *= xsi99; \
        } \
}
#define RE_DISTRIBUTE_1D(x99,new_x99,new_elem_rhs99,only99) \
{ \
        int     i99, j99; \
        double xsi99, l99[2]; \
        SVECT2D vec[2]; \
        for(i99=0;i99<2;i99++) { \
                if(i99==only99) continue; \
                DIR_LENGTH_VEC(x99[only99],new_x99[i99],x99[i99],vec); \
                for(j99=0;j99<2;j99++) \
                        l99[j99] = vec[j99].x*vec[j99].x + vec[j99].y*vec[j99].y; \
                xsi99 = sqrt(l99[0]/l99[1]); \
                new_elem_rhs99[only99*3] += (1.0-xsi99) * new_elem_rhs99[i99*3]; \
                new_elem_rhs99[only99*3+1] += (1.0-xsi99) * new_elem_rhs99[i99*3+1]; \
                new_elem_rhs99[only99*3+2] += (1.0-xsi99) * new_elem_rhs99[i99*3+2]; \
                new_elem_rhs99[i99*3] *= xsi99; \
                new_elem_rhs99[i99*3+1] *= xsi99; \
                new_elem_rhs99[i99*3+2] *= xsi99; \
        } \
}
#define RE_DISTRIBUTE_RES(x99,new_x99,new_elem_rhs99,only99) \
{ \
        int     i99, j99; \
        double xsi99, l99[2]; \
        SVECT2D vec[2]; \
        for(i99=0;i99<3;i99++) { \
                if(i99==only99) continue; \
                DIR_LENGTH_VEC(x99[only99],new_x99[i99],x99[i99],vec); \
                for(j99=0;j99<2;j99++) \
                        l99[j99] = vec[j99].x*vec[j99].x + vec[j99].y*vec[j99].y; \
                xsi99 = sqrt(l99[0]/l99[1]); \
                new_elem_rhs99[only99] += (1.0-xsi99) * new_elem_rhs99[i99]; \
                new_elem_rhs99[i99] *= xsi99; \
        } \
}

#define IS_ABOUT(a,b,c) (a<(1.+c/100.)*b && a>(1.-c/100.)*b ? 1 : 0)

#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define IS_OVER_UPPER_BOUND(a,b) (a >= b ? 1 : 0)

#define IS_BELOW_LOWER_BOUND(a,b) (a <= b ? 1 : 0)

#define IS_BTWN_BOUNDS(val,low,up) ((((low) < (val)) && ((val) < (up))) ? 1 : 0)

#define IS_OUT_BOUNDS(val,low,up) (((val) < (low) || (up) < (val)) ? 1 : 0)

#define IS_WTHN_BOUNDS(val,low,up) ((val) >= (low) ? ((val) <= (up) ? 1 : 0) : 0)

#define CELL_AVG(a) ((a[0]+a[1]+a[2])/3.0)

#define CELL_AVG_X(a) ((a[0].x+a[1].x+a[2].x)/3.0)

#define CELL_AVG_Y(a) ((a[0].y+a[1].y+a[2].y)/3.0)

#define DIST_2D(a,b) sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y))

#define DIST_3D(a,b) sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z))

#define POS_IP(result,src1,src2) \
{ \
        result = src1[0] * src2[0]; \
        result += src1[1] * src2[1]; \
}

#define POS_NORM(result,src) \
{ \
        POS_IP(result, src, src); \
        result = sqrt(result); \
}

/* macros for debugging/diagnostics  */
#define DEBUG_GENERAL (mask_debug & 15)
#define DEBUG_SOLV   ((mask_debug>>4) & 15)
#define DEBUG_FE     ((mask_debug>>8) & 15)
#define DEBUG_GRID   ((mask_debug>>12) & 15)
#define DEBUG_COMM   ((mask_debug>>16) & 15)
#define DEBUG_ELEM   ((mask_debug>>20) & 15)
#define DEBUG_NODE   ((mask_debug>>24) & 15)
#define DEBUG_OTHER  ((mask_debug>>28) & 15)
