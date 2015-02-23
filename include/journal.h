#include <stdio.h>
#include <stdlib.h>

void journal(FILE *in_f, FILE *out_f, FILE *jou_f, const char *message, char *data)
{
  const int bdim = 400;
  char buff[bdim];

  fprintf(out_f,"\n%s",message);
  fflush(out_f);
  fprintf(jou_f,"%s\n",message);
  if (in_f != stdin)
  {
    do
    {
      fgets(buff,bdim,in_f);
    } while (buff[0] == '#');
  } else
    fgets(buff,bdim,in_f);
  sscanf(buff,"%s",data);
  fprintf(jou_f,"%s\n",data);
  if (in_f != stdin) fprintf(out_f,"%s",data);
  fflush(out_f);
}

void journal(FILE *in_f, FILE *out_f, FILE *jou_f, const char *message, int &data)
{
  const int bdim = 400;
  char buff[bdim];

  fprintf(out_f,"\n%s",message);
  fflush(out_f);
  fprintf(jou_f,"%s\n",message);
  if (in_f != stdin)
  {
    do
    {
      fgets(buff,bdim,in_f);
    } while (buff[0] == '#');
  } else
    fgets(buff,bdim,in_f);
  sscanf(buff,"%d",&data);
  fprintf(jou_f,"%d\n",data);
  if (in_f != stdin) fprintf(out_f,"%d",data);
  fflush(out_f);
}

void journal(FILE *in_f, FILE *out_f, FILE *jou_f, const char *message, double &data)
{
  const int bdim = 400;
  char buff[bdim];

  fprintf(out_f,"\n%s",message);
  fflush(out_f);
  fprintf(jou_f,"%s\n",message);
  if (in_f != stdin)
  {
    do
    {
      fgets(buff,bdim,in_f);
    } while (buff[0] == '#');
  } else
    fgets(buff,bdim,in_f);
  sscanf(buff,"%lg",&data);
  fprintf(jou_f,"%lg\n",data);
  if (in_f != stdin) fprintf(out_f,"%lg",data);
  fflush(out_f);
}

void journal(FILE *in_f, FILE *out_f, FILE *jou_f, const char *message, float &data)
{
  const int bdim = 400;
  char buff[bdim];

  fprintf(out_f,"\n%s",message);
  fflush(out_f);
  fprintf(jou_f,"%s\n",message);
  if (in_f != stdin)
  {
    do
    {
      fgets(buff,bdim,in_f);
    } while (buff[0] == '#');
  } else
    fgets(buff,bdim,in_f);
  sscanf(buff,"%g",&data);
  fprintf(jou_f,"%g\n",data);
  if (in_f != stdin) fprintf(out_f,"%g",data);
  fflush(out_f);
}

void journal(FILE *in_f, FILE *out_f, FILE *jou_f, const char *message, char &data)
{
  const int bdim = 400;
  char buff[bdim];

  fprintf(out_f,"\n%s",message);
  fflush(out_f);
  fprintf(jou_f,"%s\n",message);
  if (in_f != stdin)
  {
    do
    {
      fgets(buff,bdim,in_f);
    } while (buff[0] == '#');
  } else
    fgets(buff,bdim,in_f);
  sscanf(buff,"%c",&data);
  fprintf(jou_f,"%c\n",data);
  if (in_f != stdin) fprintf(out_f,"%c",data);
  fflush(out_f);
}

