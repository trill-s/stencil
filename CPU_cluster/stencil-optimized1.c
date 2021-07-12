#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "common.h"

const char* version_name = "Optimized version";

/* your implementation */
void create_dist_grid(dist_grid_info_t *grid_info, int stencil_type) {
//   int npx,npy,npz;
   int myid,procs;
   int nx;
   int ny;
   int nz;
   int nx_global;
   int ny_global;
   int nz_global;
   int offset_x=0;
   int offset_y=0;
   int offset_z=0;
   int residux=0;
   int residuy=0;
   int residuz=0;

   myid=grid_info->p_id;
   procs=grid_info->p_num;
   int slices=cbrt(procs); /*slices on each edge*/
   nx_global=grid_info->global_size_x;
   ny_global=grid_info->global_size_y;
   nz_global=grid_info->global_size_z;

   nx=nx_global/slices;
   ny=ny_global/slices;
   nz=nz_global/slices;
   residux=nx_global%slices; 
   residuy=ny_global%slices;
   residuz=nz_global%slices;

   int x_element;
   int y_element;
   int z_element;

   z_element=myid/(slices*slices);
   y_element=(myid/slices)%slices;
   x_element=myid%slices;

   if (x_element<residux){
    nx=nx+1;
    offset_x=x_element*nx;
   }   
   else{
    offset_x=x_element*nx+residux;
   }

   if (y_element<residuy){
    ny=ny+1;
    offset_y=y_element*ny;
   }   
   else{
    offset_y=y_element*ny+residuy;
   }
   
   if (z_element<residuz){
    nz=nz+1;
    offset_z=z_element*nz;
   }   
   else{
    offset_z=z_element*nz+residuz;
   }

   /*
   if(myid<residuz) {
	nz=nz+1;
	offset_z=myid*nz;
   }
   else
   {
	offset_z=myid*nz+residuz;
   }
   */

   grid_info->local_size_x = nx;
   grid_info->local_size_y = ny;
   grid_info->local_size_z = nz;
   grid_info->offset_x = offset_x;
   grid_info->offset_y = offset_y;
   grid_info->offset_z = offset_z;
   grid_info->halo_size_x = 1;
   grid_info->halo_size_y = 1;
   grid_info->halo_size_z = 1;
   MPI_Barrier(MPI_COMM_WORLD);
	
/*    puts("not implemented");
    exit(1);*/



}

/* your implementation */
void destroy_dist_grid(dist_grid_info_t *grid_info) {

}

/* your implementation */
ptr_t stencil_7(ptr_t grid, ptr_t aux, const dist_grid_info_t *grid_info, int nt) {
    int pid = grid_info->p_id;
    int procs = grid_info->p_num;
    ptr_t buffer[2] = {grid, aux};
    int x_start = grid_info->halo_size_x, x_end = grid_info->local_size_x + grid_info->halo_size_x;
    int y_start = grid_info->halo_size_y, y_end = grid_info->local_size_y + grid_info->halo_size_y;
    int z_start = grid_info->halo_size_z, z_end = grid_info->local_size_z + grid_info->halo_size_z;
    int ldx = grid_info->local_size_x + 2 * grid_info->halo_size_x;
    int ldy = grid_info->local_size_y + 2 * grid_info->halo_size_y;
    int ldz = grid_info->local_size_z + 2 * grid_info->halo_size_z;

    int t;
    for(t = 0; t < nt; ++t) {
//bak        cptr_t a0 = buffer[t % 2];
        ptr_t a0 = buffer[t % 2];

//	printf("t=%d\n",t);
    exchange_boundary_x(a0,grid_info);
    exchange_boundary_y(a0,grid_info);
	exchange_boundary_z(a0,grid_info);


        int z, y, x;
        ptr_t a1 = buffer[(t + 1) % 2];
        for(z = z_start; z < z_end; ++z) {
            for(y = y_start; y < y_end; ++y) {
                for(x = x_start; x < x_end; ++x) {
                    a1[INDEX(x, y, z, ldx, ldy)] \
                        = ALPHA_ZZZ * a0[INDEX(x, y, z, ldx, ldy)] \
                        + ALPHA_NZZ * a0[INDEX(x-1, y, z, ldx, ldy)] \
                        + ALPHA_PZZ * a0[INDEX(x+1, y, z, ldx, ldy)] \
                        + ALPHA_ZNZ * a0[INDEX(x, y-1, z, ldx, ldy)] \
                        + ALPHA_ZPZ * a0[INDEX(x, y+1, z, ldx, ldy)] \
                        + ALPHA_ZZN * a0[INDEX(x, y, z-1, ldx, ldy)] \
                        + ALPHA_ZZP * a0[INDEX(x, y, z+1, ldx, ldy)];
                }
            }
        }
    }
    return buffer[nt % 2];
//    return grid;
}

/* your implementation */
ptr_t stencil_27(ptr_t grid, ptr_t aux, const dist_grid_info_t *grid_info, int nt) {
    ptr_t buffer[2] = {grid, aux};
    int x_start = grid_info->halo_size_x, x_end = grid_info->local_size_x + grid_info->halo_size_x;
    int y_start = grid_info->halo_size_y, y_end = grid_info->local_size_y + grid_info->halo_size_y;
    int z_start = grid_info->halo_size_z, z_end = grid_info->local_size_z + grid_info->halo_size_z;
    int ldx = grid_info->local_size_x + 2 * grid_info->halo_size_x;
    int ldy = grid_info->local_size_y + 2 * grid_info->halo_size_y;
    int ldz = grid_info->local_size_z + 2 * grid_info->halo_size_z;
    int t;
    for(t = 0; t < nt; ++t) {
//bak        cptr_t a0 = buffer[t % 2];
        ptr_t a0 = buffer[t % 2];

//	printf("t=%d\n",t);
    exchange_boundary_x(a0,grid_info);
    exchange_boundary_y(a0,grid_info);
	exchange_boundary_z(a0,grid_info);


        int z, y, x;
        ptr_t a1 = buffer[(t + 1) % 2];
        for(z = z_start; z < z_end; ++z) {
            for(y = y_start; y < y_end; ++y) {
                for(x = x_start; x < x_end; ++x) {
                    a1[INDEX(x, y, z, ldx, ldy)] \
                        = ALPHA_ZZZ * a0[INDEX(x, y, z, ldx, ldy)] \
                        + ALPHA_NZZ * a0[INDEX(x-1, y, z, ldx, ldy)] \
                        + ALPHA_PZZ * a0[INDEX(x+1, y, z, ldx, ldy)] \
                        + ALPHA_ZNZ * a0[INDEX(x, y-1, z, ldx, ldy)] \
                        + ALPHA_ZPZ * a0[INDEX(x, y+1, z, ldx, ldy)] \
                        + ALPHA_ZZN * a0[INDEX(x, y, z-1, ldx, ldy)] \
                        + ALPHA_ZZP * a0[INDEX(x, y, z+1, ldx, ldy)] \
                        + ALPHA_NNZ * a0[INDEX(x-1, y-1, z, ldx, ldy)] \
                        + ALPHA_PNZ * a0[INDEX(x+1, y-1, z, ldx, ldy)] \
                        + ALPHA_NPZ * a0[INDEX(x-1, y+1, z, ldx, ldy)] \
                        + ALPHA_PPZ * a0[INDEX(x+1, y+1, z, ldx, ldy)] \
                        + ALPHA_NZN * a0[INDEX(x-1, y, z-1, ldx, ldy)] \
                        + ALPHA_PZN * a0[INDEX(x+1, y, z-1, ldx, ldy)] \
                        + ALPHA_NZP * a0[INDEX(x-1, y, z+1, ldx, ldy)] \
                        + ALPHA_PZP * a0[INDEX(x+1, y, z+1, ldx, ldy)] \
                        + ALPHA_ZNN * a0[INDEX(x, y-1, z-1, ldx, ldy)] \
                        + ALPHA_ZPN * a0[INDEX(x, y+1, z-1, ldx, ldy)] \
                        + ALPHA_ZNP * a0[INDEX(x, y-1, z+1, ldx, ldy)] \
                        + ALPHA_ZPP * a0[INDEX(x, y+1, z+1, ldx, ldy)] \
                        + ALPHA_NNN * a0[INDEX(x-1, y-1, z-1, ldx, ldy)] \
                        + ALPHA_PNN * a0[INDEX(x+1, y-1, z-1, ldx, ldy)] \
                        + ALPHA_NPN * a0[INDEX(x-1, y+1, z-1, ldx, ldy)] \
                        + ALPHA_PPN * a0[INDEX(x+1, y+1, z-1, ldx, ldy)] \
                        + ALPHA_NNP * a0[INDEX(x-1, y-1, z+1, ldx, ldy)] \
                        + ALPHA_PNP * a0[INDEX(x+1, y-1, z+1, ldx, ldy)] \
                        + ALPHA_NPP * a0[INDEX(x-1, y+1, z+1, ldx, ldy)] \
                        + ALPHA_PPP * a0[INDEX(x+1, y+1, z+1, ldx, ldy)];
                }
            }
        }
    }
    return buffer[nt % 2];
//    return grid;
}

void exchange_boundary_x(ptr_t a0,const dist_grid_info_t *grid_info)
{
   int i,j,k,i1,nsize;
   int nx,ny,nz;
   int off_block;
   int left,right;
   int myid,procs;
   int XLAP;
   MPI_Status status;

   nx=grid_info->local_size_x;
   ny=grid_info->local_size_y;
   nz=grid_info->local_size_z;
   myid=grid_info->p_id;
   procs=grid_info->p_num;
   XLAP=grid_info->halo_size_x;

   ptr_t tmp_send1= (ptr_t)malloc(sizeof(data_t)*XLAP*ny*nz);
   ptr_t tmp_send2= (ptr_t)malloc(sizeof(data_t)*XLAP*ny*nz);
   ptr_t tmp_recv1= (ptr_t)malloc(sizeof(data_t)*XLAP*ny*nz);
   ptr_t tmp_recv2= (ptr_t)malloc(sizeof(data_t)*XLAP*ny*nz);


   int ldx = grid_info->local_size_x + 2 * grid_info->halo_size_x;
   int ldy = grid_info->local_size_y + 2 * grid_info->halo_size_y;
   int ldz = grid_info->local_size_z + 2 * grid_info->halo_size_z;
   nsize=XLAP*ny*nz;
   for(k=1;k<=nz;k++){
   for(j=1;j<=ny;j++){
    for(i=1;i<=XLAP;i++){
        i1=(k-1)*ny+j-1;
//      tmp_send1[k1]=a0[i+j*ldx+k*ldx*ldy];
//  tmp_send2[k1]=a0[i+j*ldx+(nz-ZLAP+k)*ldx*ldy];
        tmp_send1[i1]=a0[INDEX(i, j, k, ldx, ldy)];
    tmp_send2[i1]=a0[INDEX((nx-XLAP+i), j, k, ldx, ldy)];
    }
   }
   }

   if(grid_info->offset_x==0){
    left=MPI_PROC_NULL;
   }
   else
   {
    left=myid-1;
   }

   if(grid_info->offset_x==(grid_info->global_size_x-nx)){
    right=MPI_PROC_NULL;
   }
   else
   {
    right=myid+1;
   }

   MPI_Sendrecv(tmp_send1,nsize,MPI_DOUBLE,left,9000,tmp_recv2,nsize,MPI_DOUBLE,right,9000,MPI_COMM_WORLD,&status);
   MPI_Sendrecv(tmp_send2,nsize,MPI_DOUBLE,right,8000,tmp_recv1,nsize,MPI_DOUBLE,left,8000,MPI_COMM_WORLD,&status);

  if(grid_info->offset_x!=0){
   for(k=1;k<=nz;k++){
    for(j=1;j<=ny;j++){
     for(i=1;i<=XLAP;i++){
    i1=(k-1)*ny+j-1;
    a0[INDEX(i-XLAP, j, k, ldx, ldy)]=tmp_recv1[i1];
     }
    }
   }
  }


  if(grid_info->offset_x!=(grid_info->global_size_x-nx)){
   for(k=1;k<=nz;k++){
    for(j=1;j<=ny;j++){
     for(i=1;i<=XLAP;i++){
    i1=(k-1)*ny+j-1;
    a0[INDEX(nx+i, j, k, ldx, ldy)]=tmp_recv2[i1];
     }
    }
   }
  }

   free(tmp_send1);
   free(tmp_send2);
   free(tmp_recv1);
   free(tmp_recv2);
  MPI_Barrier(MPI_COMM_WORLD);
//  printf("just end MPI Comm\n");

}

void exchange_boundary_y(ptr_t a0,const dist_grid_info_t *grid_info)
{
   int i,j,k,j1,nsize;
   int nx,ny,nz;
   int off_block;
   int left,right;
   int myid,procs;
   int YLAP;
   MPI_Status status;

   nx=grid_info->local_size_x;
   ny=grid_info->local_size_y;
   nz=grid_info->local_size_z;
   myid=grid_info->p_id;
   procs=grid_info->p_num;
   YLAP=grid_info->halo_size_y;

   ptr_t tmp_send1= (ptr_t)malloc(sizeof(data_t)*nx*YLAP*nz);
   ptr_t tmp_send2= (ptr_t)malloc(sizeof(data_t)*nx*YLAP*nz);
   ptr_t tmp_recv1= (ptr_t)malloc(sizeof(data_t)*nx*YLAP*nz);
   ptr_t tmp_recv2= (ptr_t)malloc(sizeof(data_t)*nx*YLAP*nz);


   int ldx = grid_info->local_size_x + 2 * grid_info->halo_size_x;
   int ldy = grid_info->local_size_y + 2 * grid_info->halo_size_y;
   int ldz = grid_info->local_size_z + 2 * grid_info->halo_size_z;
   nsize=nx*YLAP*nz;
   for(k=1;k<=nz;k++){
   for(j=1;j<=YLAP;j++){
    for(i=1;i<=nx;i++){
        j1=(k-1)*nx+i-1;
//      tmp_send1[k1]=a0[i+j*ldx+k*ldx*ldy];
//  tmp_send2[k1]=a0[i+j*ldx+(nz-ZLAP+k)*ldx*ldy];
        tmp_send1[j1]=a0[INDEX(i, j, k, ldx, ldy)];
    tmp_send2[j1]=a0[INDEX(i, ny-YLAP+j, k, ldx, ldy)];
    }
   }
   }

   if(grid_info->offset_y==0){
    left=MPI_PROC_NULL;
   }
   else
   {
    left=myid-cbrt(procs);
   }

   if(grid_info->offset_y==grid_info->global_size_y-ny){
    right=MPI_PROC_NULL;
   }
   else
   {
    right=myid+cbrt(procs);
   }

//   printf("before sendrecv myid=%d\n",myid);

   MPI_Sendrecv(tmp_send1,nsize,MPI_DOUBLE,left,9000,tmp_recv2,nsize,MPI_DOUBLE,right,9000,MPI_COMM_WORLD,&status);
   MPI_Sendrecv(tmp_send2,nsize,MPI_DOUBLE,right,8000,tmp_recv1,nsize,MPI_DOUBLE,left,8000,MPI_COMM_WORLD,&status);

  if(grid_info->offset_y!=0){
   for(k=1;k<=nz;k++){
    for(j=1;j<=YLAP;j++){
     for(i=1;i<=nx;i++){
    j1=(k-1)*nx+i-1;
//  a0[i+j*ldx+(k-ZLAP)*ldx*ldy]=tmp_recv1[k1];
    a0[INDEX(i, (j-YLAP), k, ldx, ldy)]=tmp_recv1[j1];
     }
    }
   }
  }


  if(grid_info->offset_y!=grid_info->global_size_y-ny){
   for(k=1;k<=nz;k++){
    for(j=1;j<=YLAP;j++){
     for(i=1;i<=nx;i++){
    j1=(k-1)*nx+i-1;
//  a0[i+j*ldx+(nz+k)*ldx*ldy]=tmp_recv2[k1];
    a0[INDEX(i, (ny+j), k, ldx, ldy)]=tmp_recv2[j1];
     }
    }
   }
  }


   free(tmp_send1);
   free(tmp_send2);
   free(tmp_recv1);
   free(tmp_recv2);
  MPI_Barrier(MPI_COMM_WORLD);
//  printf("just end MPI Comm\n");

}


void exchange_boundary_z(ptr_t a0,const dist_grid_info_t *grid_info)
{
   int i,j,k,k1,nsize;
   int nx,ny,nz;
   int off_block;
   int left,right;
   int myid,procs;
   int ZLAP;
   MPI_Status status;

   nx=grid_info->local_size_x;
   ny=grid_info->local_size_y;
   nz=grid_info->local_size_z;
   myid=grid_info->p_id;
   procs=grid_info->p_num;
   ZLAP=grid_info->halo_size_z;


   ptr_t tmp_send1= (ptr_t)malloc(sizeof(data_t)*nx*ny*ZLAP);
   ptr_t tmp_send2= (ptr_t)malloc(sizeof(data_t)*nx*ny*ZLAP);
   ptr_t tmp_recv1= (ptr_t)malloc(sizeof(data_t)*nx*ny*ZLAP);
   ptr_t tmp_recv2= (ptr_t)malloc(sizeof(data_t)*nx*ny*ZLAP);


   int ldx = grid_info->local_size_x + 2 * grid_info->halo_size_x;
   int ldy = grid_info->local_size_y + 2 * grid_info->halo_size_y;
   int ldz = grid_info->local_size_z + 2 * grid_info->halo_size_z;
   nsize=nx*ny*ZLAP;
   for(k=1;k<=ZLAP;k++){
   for(j=1;j<=ny;j++){
    for(i=1;i<=nx;i++){
        k1=(k-1)*nx*ZLAP+(j-1)*nx+i-1;
//    	tmp_send1[k1]=a0[i+j*ldx+k*ldx*ldy];
//	tmp_send2[k1]=a0[i+j*ldx+(nz-ZLAP+k)*ldx*ldy];
    	tmp_send1[k1]=a0[INDEX(i, j, k, ldx, ldy)];
	tmp_send2[k1]=a0[INDEX(i, j, (nz-ZLAP+k), ldx, ldy)];
    }
   }
   }

   if(grid_info->offset_z==0){
	left=MPI_PROC_NULL;
   }
   else
   {
	left=myid-cbrt(procs)*cbrt(procs);
   }

   if(grid_info->offset_z==grid_info->global_size_z-nz){
	right=MPI_PROC_NULL;
   }
   else
   {
	right=myid+cbrt(procs)*cbrt(procs);
   }

//   printf("before sendrecv myid=%d\n",myid);

   printf("myid=%d, offset=%d, left=%d, right=%d\n",myid, grid_info->offset_x, left, right);

   MPI_Sendrecv(tmp_send1,nsize,MPI_DOUBLE,left,9000,tmp_recv2,nsize,MPI_DOUBLE,right,9000,MPI_COMM_WORLD,&status);
   MPI_Sendrecv(tmp_send2,nsize,MPI_DOUBLE,right,8000,tmp_recv1,nsize,MPI_DOUBLE,left,8000,MPI_COMM_WORLD,&status);

  if(grid_info->offset_z!=0){
   for(k=1;k<=ZLAP;k++){
    for(j=1;j<=ny;j++){
     for(i=1;i<=nx;i++){
	k1=(j-1)*nx+i-1;
//	a0[i+j*ldx+(k-ZLAP)*ldx*ldy]=tmp_recv1[k1];
	a0[INDEX(i, j, (k-ZLAP), ldx, ldy)]=tmp_recv1[k1];
     }
    }
   }
  }


  if(grid_info->offset_z!=(grid_info->global_size_z-nz)){
   for(k=1;k<=ZLAP;k++){
    for(j=1;j<=ny;j++){
     for(i=1;i<=nx;i++){
	k1=(k-1)*nx*ny+(j-1)*nx+i-1;
//	a0[i+j*ldx+(nz+k)*ldx*ldy]=tmp_recv2[k1];
	a0[INDEX(i, j, (nz+k), ldx, ldy)]=tmp_recv2[k1];
     }
    }
   }
  }
// if(myid==1){
//  printf("nx=%d,ny=%d,nz=%d,offset_z=%d\n",nx,ny,nz,grid_info->offset_z);
//  int x,y,z;
//   for(int k=grid_info->offset_z;k<=grid_info->offset_z+nz+1;k++){
//     if(k==128 || k==129 || k==256 ){
//       for(int j=1;j<=ny;j++){
//        for(int i=1;i<=nx;i++){
//          x=i;
//          y=j;
//          z=k-grid_info->offset_z;
//          printf("k,j,i,a0=%d\t,%d\t,%d\t,%lf\n",k,j,i,a0[INDEX(x,y, z, ldx, ldy)]);
//        }
//      }
//     }
//   }
// }

   free(tmp_send1);
   free(tmp_send2);
   free(tmp_recv1);
   free(tmp_recv2);
  MPI_Barrier(MPI_COMM_WORLD);
//  printf("just end MPI Comm\n");

}


/*
void exchange_corner_x(ptr_t a0,const dist_grid_info_t *grid_info)
{   
    MPI_Status status;
    int left, right;
    int myid = grid_info->p_id;
    int procs = grid_info->p_num;
    int ldx = grid_info->local_size_x + 2 * grid_info->halo_size_x;
    int ldy = grid_info->local_size_y + 2 * grid_info->halo_size_y;
    int ldz = grid_info->local_size_z + 2 * grid_info->halo_size_z;
    int slices = cbrt(procs);

    int nx=grid_info->local_size_x;
    int ny=grid_info->local_size_y;
    int nz=grid_info->local_size_z;

    float send1[4];
    send1[0]=a0[INDEX(1, 1, 1, ldx, ldy)];
    send1[1]=a0[INDEX(1, ny, 1, ldx, ldy)];
    send1[2]=a0[INDEX(1, 1, nz, ldx, ldy)];
    send1[3]=a0[INDEX(1, ny, nz, ldx, ldy)];

    float send2[4];
    send2[0]=a0[INDEX(nx, 1, 1, ldx, ldy)];
    send2[1]=a0[INDEX(nx, ny, 1, ldx, ldy)];
    send2[2]=a0[INDEX(nx, 1, nz, ldx, ldy)];
    send2[3]=a0[INDEX(nx, ny, nz, ldx, ldy)];

    float recv1[4];
    float recv2[4];


    if(grid_info->offset_x==0){
        left=MPI_PROC_NULL;
    }
    else
    {
        left=myid-1;
    }

    if(grid_info->offset_x==(grid_info->global_size_x-nx)){
        right=MPI_PROC_NULL;
    }
   else
   {
    right=myid+1;
    }

    MPI_Sendrecv(send1,4,MPI_DOUBLE,left,9000,recv2,4,MPI_DOUBLE,right,9000,MPI_COMM_WORLD,&status);
    MPI_Sendrecv(send2,4,MPI_DOUBLE,right,8000,recv1,4,MPI_DOUBLE,left,8000,MPI_COMM_WORLD,&status);

    if(right!=MPI_PROC_NULL){
        a0[nx+1,0,0,ldx,ldy]=recv2[0];
        a0[nx+1,ny+1,0,ldx,ldy]=recv2[1];
        a0[nx+1,0,nz+1,ldx,ldy]=recv2[2];
        a0[nx+1,ny+1,nz+1,ldx,ldy]=recv2[3];
    }
    

    if(left!=MPI_PROC_NULL){
        a0[0,0,0,ldx,ldy]=recv1[0];
        a0[0,ny+1,0,ldx,ldy]=recv1[1];
        a0[0,0,nz+1,ldx,ldy]=recv1[2];
        a0[0,ny+1,nz+1,ldx,ldy]=recv1[3];
    }
    MPI_Barrier(MPI_COMM_WORLD);

}
*/