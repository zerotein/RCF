#include "atoms.h"

int zout_to_prediction(Translation* zout_translations, Rotation *zout_rotations, Rotation *rotations, int number_of_angles, Prediction_Complex 
*prediction_complex, int index)
{
  int i;
  
(*prediction_complex).index_translation=(zout_translations[index]).x+(zout_translations[index]).y*hout_arguments.N+(zout_translations[index]).z*hout_arguments.N*hout_arguments.N;
  for(i=0;i<number_of_angles;i++){
    if((rotations[i]).a == (zout_rotations[index]).a &&
       (rotations[i]).b == (zout_rotations[index]).b &&
       (rotations[i]).c == (zout_rotations[index]).c){
      (*prediction_complex).index_rotation=i;
      break;
    }
  }
  return 1;
}

int read_zout_file(char* zout_file_name, Translation *translations, Rotation *rotations, double *score)
{
  FILE *zout_file;
  char line[1024],junk[1024];
  int i=1;

  hout_arguments.spacing=PREDEFINED_SPACING;
  if(zout_file=fopen(zout_file_name,"r")){
    while ((fgets(line,1024,zout_file))!=NULL){
      if(i == 1)
	sscanf(line, "%d%lf", 
	       &(hout_arguments.N), &(hout_arguments.spacing));
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
	sscanf(line, "%lf%lf%lf%d%d%d%lf", 
	       &((*(rotations+i-5)).a),
	       &((*(rotations+i-5)).b),
	       &((*(rotations+i-5)).c), 
	       &((*(translations+i-5)).x),
	       &((*(translations+i-5)).y),
	       &((*(translations+i-5)).z),
	       &(*(score+i-5)));
      }
      i++;
    }
    fclose(zout_file);
    return i-4;
  }
  else{
    fprintf(stderr, "ERROR: Can't open zout file \"%s\".\n\n", zout_file_name);
    exit (-1);
  }
}

void usage (char* name) {
  fprintf (stderr, "\nUSAGE:\n\n");
  fprintf (stderr, "   %s zout_file_name rotation_sampling_file zdock_version [hout_file_name]\n\n", name);
  fprintf (stderr, "   This program translates zdock output file to contact frequency calculation input file format.\n");
  fprintf (stderr, "   If hout_file_name isn't specified, zout_file_name.hout will be used.\n\n");
  exit(1);
}

int main (int argc, char* argv[]) 
{
  char *hout_file_name, *zout_file_name, *temp;
  FILE *hout_file;
  int  number_of_zout_lines, number_of_angles, i;
  int x, y, z;
  double *zout_score;
  Prediction_Complex *prediction_complex;
  Rotation *zout_rotations,*rotations ;
  Translation *zout_translations;

  if(argc == 4){
    zout_file_name=argv[1];
    INIT_ARRAY(hout_file_name, strlen(zout_file_name)+5);
    strcpy (hout_file_name, zout_file_name);
    hout_file_name[strlen(hout_file_name)]='.';
    hout_file_name[strlen(hout_file_name)]='h';
    hout_file_name[strlen(hout_file_name)]='o';
    hout_file_name[strlen(hout_file_name)]='u';
    hout_file_name[strlen(hout_file_name)]='t';
    hout_file_name[strlen(hout_file_name)]=0;
    temp=argv[2];
    strcpy(hout_arguments.angle_set_file_name, temp);
    temp=argv[3];
    strcpy(hout_arguments.version, temp);
  }
  else if(argc == 5){
    zout_file_name=argv[1];
    hout_file_name=argv[4];
    temp=argv[2];
    strcpy(hout_arguments.angle_set_file_name, temp);
    temp=argv[3];
    strcpy(hout_arguments.version, temp);
  }
  else{
    usage(argv[0]);
  }

  number_of_zout_lines=file_size(zout_file_name)-4;
  INIT_ARRAY(zout_translations, number_of_zout_lines);
  INIT_ARRAY(zout_rotations, number_of_zout_lines);
  INIT_ARRAY(zout_score, number_of_zout_lines);
  read_zout_file(zout_file_name, zout_translations, zout_rotations, zout_score);

  number_of_angles=file_size(hout_arguments.angle_set_file_name);
  INIT_ARRAY(rotations, number_of_angles);
  read_angle_set_file(rotations, number_of_angles, hout_arguments.angle_set_file_name);

  INIT_ARRAY(prediction_complex, 1);
  
  hout_file=fopen(hout_file_name, "w");
  fprintf(hout_file, "%d\t%3.1f\t0.0\t1\t%s\t%s\n",hout_arguments.N, hout_arguments.spacing, hout_arguments.version, 
hout_arguments.angle_set_file_name);
  fprintf(hout_file, "%f\t%f\t%f\n", (hout_arguments.random_orientation).a, (hout_arguments.random_orientation).b, 
(hout_arguments.random_orientation).c);
  fprintf(hout_file, "%s\t%.3f\t%.3f\t%.3f\n",hout_arguments.receptor_file_name, (hout_arguments.receptor_center_of_mass).x, 
(hout_arguments.receptor_center_of_mass).y, (hout_arguments.receptor_center_of_mass).z);
  fprintf(hout_file, "%s\t%.3f\t%.3f\t%.3f\n",hout_arguments.ligand_file_name, (hout_arguments.ligand_center_of_mass).x, 
(hout_arguments.ligand_center_of_mass).y, (hout_arguments.ligand_center_of_mass).z);
    
  for(i=0;i<number_of_zout_lines;i++){
    zout_to_prediction(zout_translations, zout_rotations, rotations, number_of_angles, prediction_complex, i);
    fprintf(hout_file, "%d\t%d\t%lf\t1\n", (*prediction_complex).index_rotation, (*prediction_complex).index_translation, zout_score[i]);
  }
  return 0;
}
	      
