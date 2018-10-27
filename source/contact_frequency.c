#include "atoms.h"

int  calculate_contact_frequency(Coordinates_Complex*, Coordinates_Complex*, unsigned int*, unsigned int*, double);
int  readLISTfile(char*, int*);

int readLISTfile(char* list_file_name, int* list_of_predictions)
{
  FILE *list_file;
  char line[200],junk[200];
  int i;
  
  i=0;
  if(list_file=fopen(list_file_name,"r")){
    while ((fgets(line,200,list_file))!=NULL){
      sscanf(line, "%d", list_of_predictions+i);
      i++;
    }
    fclose (list_file);
  }
  else{
    fprintf (stderr, "Can't open list file %s\n\n", list_file_name);
    exit ( 0 );
  }
  return i;
}

int calculate_contact_frequency(Coordinates_Complex* receptor_coordinates, Coordinates_Complex* ligand_coordinates, unsigned int* frequency_receptor, unsigned int* frequency_ligand, double radius_cutoff)
{
  int i, j;

  for (i=0;i<(*receptor_coordinates).number_of_atoms;i++){
    for (j=0;j<(*ligand_coordinates).number_of_atoms;j++){
      if(pow((((*receptor_coordinates).coordinates[i]).x-((*ligand_coordinates).coordinates[j]).x),2)+
	 pow((((*receptor_coordinates).coordinates[i]).y-((*ligand_coordinates).coordinates[j]).y),2)+
	 pow((((*receptor_coordinates).coordinates[i]).z-((*ligand_coordinates).coordinates[j]).z),2)<=radius_cutoff*radius_cutoff){
	(*(frequency_receptor+i))++;
	(*(frequency_ligand+j))++;
      }
    }
  }

  return i;
}

void usage (char* name) {
  fprintf (stderr, "\nUSEAGE:\n\n");
  fprintf (stderr, "   %s -i hout_file [-o output_nametag]\n\n", name);
  fprintf (stderr, "   This program calculates the contact frequency of all atoms using ZDOCK results.\n");
  fprintf (stderr, "   Example: %s -i 1a2k.hout -o RCF\n\n", name);
  exit(1);
}

int main (int argc, char* argv[]) 
{
  char *hout_file_name, *receptor_file_name, *ligand_file_name, *output_nametag, *output_receptor_file_name, *output_ligand_file_name, 
*list_file_name, temp_save[200], temp_save2[200];
  FILE *output_receptor_file, *output_ligand_file, *receptor_file, *ligand_file;
  int default_L=1, default_x=1, default_l=1, default_r=1, default_i=1, default_o=1, default_c=1;
  int i,j;
  int number_of_angles, number_of_list_lines, number_of_hout_lines, size_of_receptor, size_of_ligand, number_of_predictions, surface;
  int *list_of_predictions;
  unsigned int *contact_frequency_receptor, *contact_frequency_ligand;
  double radius_cutoff, contact_frequency_max, contact_frequency_min, temperature;
  Prediction_Complex  *prediction_complex;
  Position_Complex    *position_complex;
  Rotation            *rotations;
  Coordinates_Complex *ligand_coordinates_original, *ligand_coordinates, *receptor_coordinates;

  if (argc < 3)
    usage (argv[0]);

  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      switch (argv[i][1]) {
      case 'i':
	i++;
	hout_file_name=argv[i];
	default_i=0;
	break;
      case 'o':
	i++;
      	output_nametag=argv[i];
	default_o=0;
	break;
      case 'l':
	i++;
	ligand_file_name=argv[i];
	default_l=0;
	break;
      case 'r':
	i++;
	receptor_file_name=argv[i];
	default_r=0;
	break;
      case 'c':
	i++;
	radius_cutoff=atof(argv[i]);
	default_c=0;
	break;
      case 'L':
	i++;
      	list_file_name=argv[i];
	default_L=0;
	break;
      default:
	usage (argv[0]);
      }
    }
  }

  if(default_i){
    usage (argv[0]);
  }
  else{
    number_of_hout_lines=file_size(hout_file_name)-4;
    INIT_ARRAY(prediction_complex, number_of_hout_lines);
    read_hout_file(hout_file_name, prediction_complex);
  }

  if(default_l){
    get_PDB_filenames(hout_file_name);
    INIT_ARRAY(ligand_file_name, 1024);
    strcpy(ligand_file_name, hout_arguments.ligand_file_name);
  }

  if(default_r){
    get_PDB_filenames(hout_file_name);
    INIT_ARRAY(receptor_file_name, 1024);
    strcpy(receptor_file_name, hout_arguments.receptor_file_name);
  }

  else{
    if(54000>number_of_hout_lines+1){
      fprintf(stderr, "ERROR: xto exceeds limit. \n\n");
      exit (0);
    }
  } 

  if(default_c){
    radius_cutoff=6;
  }

  INIT_ARRAY(output_receptor_file_name, 2048);
  INIT_ARRAY(output_ligand_file_name, 2048);
  if(default_o){
    sprintf(output_receptor_file_name, "%s.contact_frequency", receptor_file_name);
    sprintf(output_ligand_file_name, "%s.contact_frequency", ligand_file_name);
  }
  else{
    sprintf(output_receptor_file_name, "%s.%s", receptor_file_name, output_nametag);
    sprintf(output_ligand_file_name, "%s.%s", ligand_file_name, output_nametag);
  }
  output_receptor_file=fopen(output_receptor_file_name,"w");
  output_ligand_file=fopen(output_ligand_file_name,"w");

  INIT_ARRAY(ligand_coordinates_original, 1);
  INIT_ARRAY(ligand_coordinates, 1);
  INIT_ARRAY(receptor_coordinates, 1);
  INIT_ARRAY(position_complex, 1);

  read_PDB_file(ligand_file_name, ligand_coordinates_original);
  read_PDB_file(receptor_file_name, receptor_coordinates);

  size_of_receptor=file_size(receptor_file_name);
  INIT_ARRAY(contact_frequency_receptor, size_of_receptor);
  for (i=0; i<size_of_receptor; i++){
    *(contact_frequency_receptor+i)=0;
  }
  size_of_ligand=file_size(ligand_file_name);
  INIT_ARRAY(contact_frequency_ligand, size_of_ligand);
  for (i=0; i<size_of_ligand; i++){
    *(contact_frequency_ligand+i)=0;
  }

  number_of_angles=file_size(hout_arguments.angle_set_file_name);
  INIT_ARRAY(rotations, number_of_angles);
  read_angle_set_file(rotations, number_of_angles, hout_arguments.angle_set_file_name);

  number_of_predictions=0;
  if(!default_L){
    number_of_list_lines=file_size(list_file_name);
    INIT_ARRAY(list_of_predictions, number_of_list_lines);
    readLISTfile(list_file_name, list_of_predictions);

    for (i=0;i<number_of_list_lines;i++){
      j=*(list_of_predictions+i);
      if(j>=1 && j<=2000){
	pick_prediction(position_complex, rotations, prediction_complex, j-1);
	create_atoms(ligand_coordinates_original, ligand_coordinates, position_complex);
	calculate_contact_frequency(receptor_coordinates, ligand_coordinates, contact_frequency_receptor, contact_frequency_ligand, radius_cutoff);
	number_of_predictions++;
      }
    } 
  }
  else{
    for (j=1;j<=2000;j++){
      pick_prediction(position_complex, rotations,  prediction_complex, j-1);
      create_atoms(ligand_coordinates_original, ligand_coordinates, position_complex);
      calculate_contact_frequency(receptor_coordinates, ligand_coordinates, contact_frequency_receptor, contact_frequency_ligand, radius_cutoff);
      number_of_predictions++;
    }
  }

  contact_frequency_max=*(contact_frequency_receptor);
  contact_frequency_min=*(contact_frequency_receptor);
  for (i=0;i<size_of_receptor;i++){
    if(*(contact_frequency_receptor+i)<contact_frequency_min){contact_frequency_min=*(contact_frequency_receptor+i);}
    if(*(contact_frequency_receptor+i)>contact_frequency_max){contact_frequency_max=*(contact_frequency_receptor+i);}
  }
  
  receptor_file = fopen(receptor_file_name, "r");
  for (i=0;i<size_of_receptor;i++){
    fgets(temp_save,61,receptor_file);
    fgets(temp_save2,200,receptor_file);
    temperature=((*(contact_frequency_receptor+i))-contact_frequency_min)/(contact_frequency_max-contact_frequency_min);
    if(temperature<10){
      fprintf(output_receptor_file, "%s  %1.2f\t%d\t%f\t%c\n",temp_save, temperature, *(contact_frequency_receptor+i), (*(contact_frequency_receptor+i))/(double)number_of_predictions, temp_save2[2]);
    }
    else{
      fprintf(output_receptor_file, "%s %2.2f\t%d\t%f\t%c\n",temp_save, temperature, *(contact_frequency_receptor+i), (*(contact_frequency_receptor+i))/(double)number_of_predictions, temp_save2[2]);
    }
  }
    
  contact_frequency_max=*(contact_frequency_ligand);
  contact_frequency_min=*(contact_frequency_ligand);
  for (i=0;i<size_of_ligand;i++){
    if(*(contact_frequency_ligand+i)<contact_frequency_min){contact_frequency_min=*(contact_frequency_ligand+i);}
    if(*(contact_frequency_ligand+i)>contact_frequency_max){contact_frequency_max=*(contact_frequency_ligand+i);}
  }
  
  ligand_file = fopen(ligand_file_name, "r");
  for (i=0;i<size_of_ligand;i++){
    fgets(temp_save,61,ligand_file);
    fgets(temp_save2,200,ligand_file);
    temperature=((*(contact_frequency_ligand+i))-contact_frequency_min)/(contact_frequency_max-contact_frequency_min);
    if(temperature<10){
      fprintf(output_ligand_file, "%s  %1.2f\t%d\t%f\t%c\n",temp_save, temperature, *(contact_frequency_ligand+i), (*(contact_frequency_ligand+i))/(double)number_of_predictions, temp_save2[2]);
    }
    else{
      fprintf(output_ligand_file, "%s %2.2f\t%d\t%f\t%c\n",temp_save, temperature, *(contact_frequency_ligand+i), (*(contact_frequency_ligand+i))/(double)number_of_predictions, temp_save2[2]);
    }
  }
    
  return 0;
}

