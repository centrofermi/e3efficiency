#include <TGrid.h>

void DoMerge(char* nomelista="lista",char* output="allBOLO-01.root"){
// system("ls BOLO-04-BOLO-03-2015*.root >lista");

  FILE *f = fopen(nomelista,"r");

  TFileMerger m(kFALSE);
  m.OutputFile(output);

  Int_t i=0;
  char nome[100];
  while (fscanf(f,"%s",nome)==1) {
    m.AddFile(nome);
    i++;
  }
  if (i)
    m.Merge();
}

