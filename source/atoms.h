#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>

#define INIT_ARRAY(array,size)                         \
        (*((void**)(&(array)))=(void*)malloc(sizeof(*(array))*(size)))

#define PREDEFINED_SPACING    1.2f
#define MAX_NUMBER_OF_ATOMS  32767

typedef struct Translation
{ 
  int x;
  int y;
  int z;
  int index;
} Translation;

typedef struct  Rotation
{
  double a;
  double b;
  double c;
  int    index;
}Rotation;

typedef	struct	Position_Complex
{
  double a;
  double b;
  double c;
  double x;
  double y;
  double z;
}Position_Complex;

typedef struct  Coordinates
{
  double x;
  double y;
  double z;
  int atom_index;
  int residue_index;
}Coordinates;

typedef struct  Coordinates_Complex
{
  Coordinates      coordinates[MAX_NUMBER_OF_ATOMS];
  int              number_of_atoms;
  int              number_of_residues;
}Coordinates_Complex;

typedef	struct	Prediction_Complex
{
  int    index_rotation;
  int    index_translation;
  double score;
  int    rank_within_rotation;
  double rmsd;
  char   indicator[2];
}Prediction_Complex;

typedef struct  hout_Arguments
{
  int              N;
  double           spacing;
  double           clustering_radius;
  int              max_number_of_clusters;
  char             version[4];
  char             angle_set_file_name[1024];
  Rotation         random_orientation;
  char             receptor_file_name[1024];
  Coordinates      receptor_center_of_mass;
  char             ligand_file_name[1024];
  Coordinates      ligand_center_of_mass;
}hout_Arguments;

//Override functions//
void copy_coordinates(Coordinates*, Coordinates*);
void increment(int*, int, char*);
int read_angle_set_file(Rotation*, int, char*);
int file_size(char*);
int sort_int(int*, int*, int);
int sort_double(double*, int*, int);
int sort_float(float*, int*, int);
int sort_float_partial(float*, short unsigned int*, int, int);
int get_PDB_filenames(char*);
int read_hout_file(char*, Prediction_Complex*);
int read_PDB_file(char*, Coordinates_Complex*);
int rotate_atom (Coordinates*, Coordinates*, Rotation*);
int create_atoms(Coordinates_Complex*, Coordinates_Complex*, Position_Complex*);
int pick_prediction(Position_Complex*, Rotation*, Prediction_Complex*, int);
int find_interface_by_atom(Coordinates_Complex*, Coordinates_Complex*, Coordinates_Complex*, double);
int find_interface_by_residue(Coordinates_Complex*, Coordinates_Complex*, Coordinates_Complex*, double);
int combine_coordinates(Coordinates_Complex*, Coordinates_Complex*, Coordinates_Complex*);
char* add_back_slash(char* string);

hout_Arguments hout_arguments;

void copy_coordinates(Coordinates* target, Coordinates* source)
{
  (*target).x=(*source).x;
  (*target).y=(*source).y;
  (*target).z=(*source).z;
  (*target).atom_index=(*source).atom_index;
  (*target).residue_index=(*source).residue_index;
}

void increment(int* index, int roof, char* error_message)
{
  (*index)+=1;
  if( (*index)>=roof){
    fprintf(stderr, "\nERROR: %s!!\n\n", error_message);
    exit (1);
  }
}

int file_size(char *file_name)
{
  FILE          *file;
  int           lines=0;
  char          line[1024];
  char          temp_name[1024];

  lines=0;
  if(!(file=fopen(file_name,"rt"))){
    fprintf(stderr,"\nWARNING: Can't find file \"%s\" in current directory.\n", file_name);
    sprintf(temp_name, "./%s", file_name);
    if(!(file=fopen(temp_name, "rt"))){
      fprintf(stderr,"\nERROR: Unable to open angle set file \"%s\".\n\n", file_name);
      exit(1);
    }
    strcpy(file_name, temp_name);
    fprintf(stderr,"WARNING: Automatically Switch to \"%s\"\n\n", temp_name);
  }

  while(!feof(file)){
    line[0]=0;
    if(fgets(line,1000,file)==NULL)
      break;
    lines++;
  }
  return lines;
}

int  sort_int(int* array, int* index, int size)
{
  int i, j, head, tail;

  if(size<=1){
    index[0]=0;
    return 1;
  }
  if(array[0]>array[1]){
    index[0]=0; index[1]=1;  }
  else{
    index[0]=1; index[1]=0;  }
  for (i=2; i<size; i++){
    head=0; tail=i-1; 
    while (tail-head>1){
      if(array[i]>array[index[head+(int)((tail-head)*0.618)]]){
	tail=head+(int)((tail-head)*0.618);
      }
      else{
	head=head+(int)((tail-head)*0.618);
      }
    }

    if(array[i]>array[index[head]]){
      for (j=i-1; j>=head; j--){
	index[j+1]=index[j];
      }
      index[head]=i;
    }
    else if(array[i]>array[index[tail]]){
      for (j=i-1; j>=tail; j--){
	index[j+1]=index[j];
      }
      index[tail]=i;
    }
    else{
      index[tail+1]=i;
    }
  }
  return 1;
}
int  sort_double(double* array, int* index, int size)
{
  int i, j, head, tail;

  if(size<=1){
    index[0]=0;
    return 1;
  }
  if(array[0]>array[1]){
    index[0]=0; index[1]=1;  }
  else{
    index[0]=1; index[1]=0;  }
  for (i=2; i<size; i++){
    head=0; tail=i-1;
    while (tail-head>1){
      if(array[i]>array[index[head+(int)((tail-head)*0.618)]]){
	tail=head+(int)((tail-head)*0.618);
      }
      else{
	head=head+(int)((tail-head)*0.618);
      }
    }

    if(array[i]>array[index[head]]){
      for (j=i-1; j>=head; j--){
	index[j+1]=index[j];
      }
      index[head]=i;
    }
    else if(array[i]>array[index[tail]]){
      for (j=i-1; j>=tail; j--){
	index[j+1]=index[j];
      }
      index[tail]=i;
    }
    else{
      index[tail+1]=i;
    }
  }
  return 1;
}
int  sort_float(float* array, int* index, int size)
{
  int i, j, head, tail;

  if(size<=1){
    index[0]=0;
    return 1;
  }
  if(array[0]>array[1]){
    index[0]=0; index[1]=1;  }
  else{
    index[0]=1; index[1]=0;  }
  for (i=2; i<size; i++){
    head=0; tail=i-1;  
    while (tail-head>1){
      if(array[i]>array[index[head+(int)((tail-head)*0.618)]]){
	tail=head+(int)((tail-head)*0.618);
      }
      else{
	head=head+(int)((tail-head)*0.618);
      }
    }

    if(array[i]>array[index[head]]){
      for (j=i-1; j>=head; j--){
	index[j+1]=index[j];
      }
      index[head]=i;
    }
    else if(array[i]>array[index[tail]]){
      for (j=i-1; j>=tail; j--){
	index[j+1]=index[j];
      }
      index[tail]=i;
    }
    else{
      index[tail+1]=i;
    }
  }
  return 1;
}

int  sort_float_partial(float* array, short unsigned int* index, int size, int n)
{
  int i, j, head, tail;

  if(size<=1 || n<=1){
    index[0]=0;
    return 1;
  }
  if(array[0]>array[1]){
    index[0]=0; index[1]=1;  }
  else{
    index[0]=1; index[1]=0;  }
  for (i=2; i<size; i++){

    head=0; tail=i-1;
    if(i>=n)tail=n-2;    
    while (tail-head>1){
      if(array[i]>array[index[head+(int)((tail-head)*0.618)]]){
	tail=head+(int)((tail-head)*0.618);
      }
      else{
	head=head+(int)((tail-head)*0.618);
      }
    }

    if(array[i]>array[index[head]]){
      for (j=i-1; j>=head; j--){
	index[j+1]=index[j];
      }
      index[head]=i;
    }
    else if(array[i]>array[index[tail]]){
      for (j=i-1; j>=tail; j--){
	index[j+1]=index[j];
      }
      index[tail]=i;
    }
    else{
      index[tail+1]=i;
    }
  }
  return 1;
}

int read_angle_set_file(Rotation *rotations, int number_of_rotations, char *angle_set_file_name)
{
  FILE          *angle_set_file;
  int           number_of_angles;
  char          line[1024];
  
  if(!(angle_set_file=fopen(angle_set_file_name,"rt"))){
    fprintf(stderr, "\nERROR: Can't open angle set file!!\n\n");
    exit (-1);
  }

  while(!feof(angle_set_file)){
    line[0]=0;
    if(fgets(line,1024,angle_set_file)==NULL)break;
    if(sscanf(line, "%d%lf%lf%lf", &((*(rotations)).index), &((*(rotations)).a), &((*(rotations)).b), &((*(rotations)).c))  == 4){
      rotations++;
    }
    else{
      fprintf(stderr, "\nERROR: Angle set file doesn't follow the right format - index[\\t]angle1[\\t]angle2[\\t]angle3\n\n");
      exit (-1);
    }
  }
  return 1;
}

int get_PDB_filenames(char* hout_file_name)
{
  FILE *hout_file;
  char line[1024];
  double undef;
  int i=1;
  
  if(hout_file=fopen(hout_file_name,"r")){
    while ((fgets(line,1024,hout_file))!=NULL){
      if(i == 3){
	sscanf(line, "%s%lf%lf%lf", hout_arguments.receptor_file_name, &undef, &undef, &undef);
      }
      else if(i == 4){
	sscanf(line, "%s%lf%lf%lf", hout_arguments.ligand_file_name, &undef, &undef, &undef);
	break;
      }
      i++;
    }
  }

  fclose(hout_file);
  return 1;
}

int read_hout_file(char* hout_file_name, Prediction_Complex* prediction_complex)
{
  FILE *hout_file;
  char line[1024],junk[1024];
  int i=1, t1, t2, t3;

  hout_arguments.spacing=PREDEFINED_SPACING;
  if(hout_file=fopen(hout_file_name,"r")){
    while ((fgets(line,1024,hout_file))!=NULL){
      if(i == 1)
	sscanf(line, "%d%lf%lf%d%s%s", 
	       &(hout_arguments.N), &(hout_arguments.spacing), &(hout_arguments.clustering_radius), &(hout_arguments.max_number_of_clusters), hout_arguments.version, hout_arguments.angle_set_file_name);
      else if(i == 2)
	sscanf(line, "%lf%lf%lf", 
	       &((hout_arguments.random_orientation).a), &((hout_arguments.random_orientation).b), &((hout_arguments.random_orientation).c));
      else if(i == 3)
	sscanf(line, "%s%lf%lf%lf", 
	       hout_arguments.receptor_file_name, &((hout_arguments.receptor_center_of_mass).x), &((hout_arguments.receptor_center_of_mass).y), &((hout_arguments.receptor_center_of_mass).z));
      else if(i == 4)
	sscanf(line, "%s%lf%lf%lf", 
	       hout_arguments.ligand_file_name, &((hout_arguments.ligand_center_of_mass).x), &((hout_arguments.ligand_center_of_mass).y), &((hout_arguments.ligand_center_of_mass).z));
      else{
	sscanf(line, "%d%d%lf%d%lf%s", 
	       &((*(prediction_complex+i-5)).index_rotation),
	       &((*(prediction_complex+i-5)).index_translation),
	       &((*(prediction_complex+i-5)).score), 
	       &((*(prediction_complex+i-5)).rank_within_rotation),
	       &((*(prediction_complex+i-5)).rmsd),
	       &((*(prediction_complex+i-5)).indicator));
      }
      i++;
    }
    fclose(hout_file);
    return i-4;
  }
  else{
    fprintf(stderr, "ERROR: Can't open hout file \"%s\".\n\n", hout_file_name);
    exit (-1);
  }
}
int read_PDB_file(char *pdb_file_name, Coordinates_Complex *coordinates_complex)
{
  FILE *pdb_file;
  char getline[1024], temp_name[1024];
  int current_residue, last_residue;

  if(!(pdb_file=fopen(pdb_file_name,"rt"))){
    fprintf(stderr,"\nWARNING: Can't find file \"%s\" in current directory.\n", pdb_file_name);
    sprintf(temp_name, "./%s", pdb_file_name);
    if(!(pdb_file=fopen(temp_name, "rt"))){
      fprintf(stderr,"\nERROR: Unable to open PDB file \"%s\".\n\n", pdb_file_name);
      exit(1);
    }
    strcpy(pdb_file_name, temp_name);
    fprintf(stderr,"WARNING: Automatically Switch to \"%s\"\n\n", temp_name);
  }

  last_residue=0;
 (*coordinates_complex).number_of_atoms=0;
 (*coordinates_complex).number_of_residues=0;
 while ((fgets(getline,8,pdb_file))!=NULL){
    fgets(getline,5,pdb_file);
    ((*coordinates_complex).coordinates[(*coordinates_complex).number_of_atoms]).atom_index=atoi(getline);
    fgets(getline,12,pdb_file);
    fgets(getline,5,pdb_file);
    current_residue=atoi(getline);
    if(current_residue!=last_residue){
      (*coordinates_complex).number_of_residues++;
    }
    ((*coordinates_complex).coordinates[(*coordinates_complex).number_of_atoms]).residue_index=(*coordinates_complex).number_of_residues-1;
    last_residue=current_residue;
    fgets(getline,5,pdb_file);
    fgets(getline,9,pdb_file);
    ((*coordinates_complex).coordinates[(*coordinates_complex).number_of_atoms]).x=atof(getline);
    fgets(getline,9,pdb_file);
    ((*coordinates_complex).coordinates[(*coordinates_complex).number_of_atoms]).y=atof(getline);
    fgets(getline,9,pdb_file);
    ((*coordinates_complex).coordinates[(*coordinates_complex).number_of_atoms]).z=atof(getline);

    fgets(getline,40,pdb_file);
    increment(&((*coordinates_complex).number_of_atoms), MAX_NUMBER_OF_ATOMS, "Atoms exceeds limit");
  }
  fclose(pdb_file);

  return 1;
}

int pick_prediction(Position_Complex* position_complex, Rotation *rotations, Prediction_Complex *prediction_complex, int index)
{
  int x, y, z;

  (*position_complex).a=rotations[(prediction_complex[index]).index_rotation].a;
  (*position_complex).b=rotations[(prediction_complex[index]).index_rotation].b;
  (*position_complex).c=rotations[(prediction_complex[index]).index_rotation].c;
  x=(int)(((prediction_complex[index]).index_translation)%(hout_arguments.N));
  z=(int)(((prediction_complex[index]).index_translation)/((hout_arguments.N)*(hout_arguments.N)));
  y=(int)((((prediction_complex[index]).index_translation)-x-z*((hout_arguments.N)*(hout_arguments.N)))/(hout_arguments.N));

  if (x >= hout_arguments.N/2) x -= hout_arguments.N;
  if (y >= hout_arguments.N/2) y -= hout_arguments.N;  
  if (z >= hout_arguments.N/2) z -= hout_arguments.N;  

  (*position_complex).x=(double)x*hout_arguments.spacing; 
  (*position_complex).y=(double)y*hout_arguments.spacing; 
  (*position_complex).z=(double)z*hout_arguments.spacing;

  return 1;
}

int create_atoms(Coordinates_Complex* coordinates_complex, Coordinates_Complex* coordinates_repositioned, Position_Complex* position_complex)
{
  int          i;
  Coordinates  c1, c2, c3;
  Rotation     temp_rot, temp_rot2;

  temp_rot.a=(hout_arguments.random_orientation).a;
  temp_rot.b=(hout_arguments.random_orientation).b;
  temp_rot.c=(hout_arguments.random_orientation).c;
  temp_rot2.a=(*position_complex).a;
  temp_rot2.b=(*position_complex).b;
  temp_rot2.c=(*position_complex).c;

  for (i=0;i<(*coordinates_complex).number_of_atoms;i++){
    c1.x=((*coordinates_complex).coordinates[i]).x-(hout_arguments.ligand_center_of_mass).x;
    c1.y=((*coordinates_complex).coordinates[i]).y-(hout_arguments.ligand_center_of_mass).y;
    c1.z=((*coordinates_complex).coordinates[i]).z-(hout_arguments.ligand_center_of_mass).z;

    rotate_atom(&c1, &c2, &temp_rot);
    rotate_atom(&c2, &c3, &temp_rot2);

    ((*coordinates_repositioned).coordinates[i]).x = c3.x-(*position_complex).x+(hout_arguments.receptor_center_of_mass).x;
    ((*coordinates_repositioned).coordinates[i]).y = c3.y-(*position_complex).y+(hout_arguments.receptor_center_of_mass).y;
    ((*coordinates_repositioned).coordinates[i]).z = c3.z-(*position_complex).z+(hout_arguments.receptor_center_of_mass).z;


    ((*coordinates_repositioned).coordinates[i]).atom_index = ((*coordinates_complex).coordinates[i]).atom_index;
    ((*coordinates_repositioned).coordinates[i]).residue_index = ((*coordinates_complex).coordinates[i]).residue_index;
  }
  (*coordinates_repositioned).number_of_atoms=(*coordinates_complex).number_of_atoms;
  (*coordinates_repositioned).number_of_residues=(*coordinates_complex).number_of_residues;
  return 1;
}

int rotate_atom (Coordinates *original_coordinates, Coordinates *new_coordinates, Rotation *rotation)
{
  double r11, r21, r31, r12, r22, r32, r13, r23, r33;
  double psi, phi, theta;
  psi=(*rotation).a;
  theta=(*rotation).b;
  phi=(*rotation).c;

  r11 = cos(psi)*cos(phi)  -  sin(psi)*cos(theta)*sin(phi);
  r21 = sin(psi)*cos(phi)  +  cos(psi)*cos(theta)*sin(phi);
  r31 = sin(theta)*sin(phi);

  r12 = -cos(psi)*sin(phi)  -  sin(psi)*cos(theta)*cos(phi);
  r22 = -sin(psi)*sin(phi)  +  cos(psi)*cos(theta)*cos(phi);
  r32 = sin(theta)*cos(phi);

  r13 = sin(psi)*sin(theta);
  r23 = -cos(psi)*sin(theta);
  r33 = cos(theta);

  (*new_coordinates).x = r11 * (*original_coordinates).x + r12 * (*original_coordinates).y + r13 * (*original_coordinates).z;
  (*new_coordinates).y = r21 * (*original_coordinates).x + r22 * (*original_coordinates).y + r23 * (*original_coordinates).z;
  (*new_coordinates).z = r31 * (*original_coordinates).x + r32 * (*original_coordinates).y + r33 * (*original_coordinates).z;

  return 1;
}

int find_interface_by_atom(Coordinates_Complex* moleculeA, Coordinates_Complex* moleculeB, Coordinates_Complex* interface, double interface_cutoff)
{
  int    i, j, k;
  short unsigned int    *on_interface_moleculeA, *on_interface_moleculeB;

  INIT_ARRAY(on_interface_moleculeA, (*moleculeA).number_of_atoms);
  INIT_ARRAY(on_interface_moleculeB, (*moleculeB).number_of_atoms);
  for (i=0; i<(*moleculeA).number_of_atoms; i++) on_interface_moleculeA[i]=0;
  for (i=0; i<(*moleculeB).number_of_atoms; i++) on_interface_moleculeB[i]=0;

  for (i=0; i<(*moleculeA).number_of_atoms; i++){
    for (j=0; j<(*moleculeB).number_of_atoms; j++){
      if(on_interface_moleculeA[i] && on_interface_moleculeB[j]) continue;
      if( pow((((*moleculeA).coordinates[i]).x-((*moleculeB).coordinates[j]).x),2)+
	  pow((((*moleculeA).coordinates[i]).y-((*moleculeB).coordinates[j]).y),2)+
	  pow((((*moleculeA).coordinates[i]).z-((*moleculeB).coordinates[j]).z),2)
	  <= interface_cutoff*interface_cutoff ){
	  on_interface_moleculeA[i]=1;
	  on_interface_moleculeB[j]=1;
      }
    }
  }
  k=0;
  for (i=0; i<(*moleculeA).number_of_atoms; i++){
    if(on_interface_moleculeA[i]){
      copy_coordinates( ((*interface).coordinates)+k, ((*moleculeA).coordinates)+i);
     increment(&k, MAX_NUMBER_OF_ATOMS, "Interface atoms exceeds limit");
    }
  }
  for (i=0; i<(*moleculeB).number_of_atoms; i++){
    if(on_interface_moleculeB[i]){
      copy_coordinates( ((*interface).coordinates)+k, ((*moleculeA).coordinates)+i);
      increment(&k, MAX_NUMBER_OF_ATOMS, "Interface atoms exceeds limit");
    }
  }
  (*interface).number_of_atoms=k;

  free(on_interface_moleculeA);
  free(on_interface_moleculeB);
  return 1;
}

int find_interface_by_residue(Coordinates_Complex* moleculeA, Coordinates_Complex* moleculeB, Coordinates_Complex* interface, double interface_cutoff)
{
  int    i, j, k;
  short unsigned int    *on_interface_moleculeA, *on_interface_moleculeB;

  INIT_ARRAY(on_interface_moleculeA, (*moleculeA).number_of_residues);
  INIT_ARRAY(on_interface_moleculeB, (*moleculeB).number_of_residues);
  for (i=0; i<(*moleculeA).number_of_residues; i++) on_interface_moleculeA[i]=0;
  for (i=0; i<(*moleculeB).number_of_residues; i++) on_interface_moleculeB[i]=0;

  for (i=0; i<(*moleculeA).number_of_atoms; i++){
    for (j=0; j<(*moleculeB).number_of_atoms; j++){
      if(on_interface_moleculeA[((*moleculeA).coordinates[i]).residue_index] && on_interface_moleculeB[((*moleculeB).coordinates[j]).residue_index]) continue;
      if( pow((((*moleculeA).coordinates[i]).x-((*moleculeB).coordinates[j]).x),2)+
	  pow((((*moleculeA).coordinates[i]).y-((*moleculeB).coordinates[j]).y),2)+
	  pow((((*moleculeA).coordinates[i]).z-((*moleculeB).coordinates[j]).z),2)
	  <= interface_cutoff*interface_cutoff ){
	  on_interface_moleculeA[((*moleculeA).coordinates[i]).residue_index]=1;
	  on_interface_moleculeB[((*moleculeB).coordinates[j]).residue_index]=1;
      }
    }
  }
  k=0;
  for (i=0; i<(*moleculeA).number_of_atoms; i++){
    if(on_interface_moleculeA[((*moleculeA).coordinates[i]).residue_index]){
      copy_coordinates(((*interface).coordinates)+k, ((*moleculeA).coordinates)+i);
      increment(&k, MAX_NUMBER_OF_ATOMS, "Interface atoms exceeds limit");
    }
  }
  for (i=0; i<(*moleculeB).number_of_atoms; i++){
    if(on_interface_moleculeB[((*moleculeB).coordinates[i]).residue_index]){
      copy_coordinates(((*interface).coordinates)+k, ((*moleculeB).coordinates)+i);
      increment(&k, MAX_NUMBER_OF_ATOMS, "Interface atoms exceeds limit"); 
    }
  }
  (*interface).number_of_atoms=k;

  free(on_interface_moleculeA);
  free(on_interface_moleculeB);
  return 1;
}
int combine_coordinates(Coordinates_Complex* moleculeA, Coordinates_Complex* moleculeB, Coordinates_Complex* moleculeC)
{
  int    i, k;
  k=0;
  for (i=0; i<(*moleculeA).number_of_atoms; i++){
    copy_coordinates(((*moleculeC).coordinates)+k, ((*moleculeA).coordinates)+i);
    increment(&k, MAX_NUMBER_OF_ATOMS, "Combined atoms exceeds limit");
  }
  for (i=0; i<(*moleculeB).number_of_atoms; i++){
    copy_coordinates(((*moleculeC).coordinates)+k, ((*moleculeB).coordinates)+i);
    increment(&k, MAX_NUMBER_OF_ATOMS, "Combined atoms exceeds limit");
  }
  (*moleculeC).number_of_atoms=k;
  return 1;
}

char* add_back_slash(char* string){
  int i;
  i=strlen(string);
  if(*(string+i-1) != '/')
    strcat(string,"/");
  return string;
}
