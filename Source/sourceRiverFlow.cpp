/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments,          *
 *  river networks, lakes and landscape elements.                                       *
 *  Preprocessing for traversing and defining a hierarchy of watercourse elements.      *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "classRiverFlow.h"

void GetRowCol(int elementNo, int nCols, int &row, int &col); 
void ReadFlowDirection(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
int ControlFlowDirection(int, int);
void FindUpLandFlow(DistributedElement * const Dew, int nRows, int nCols, int noData);
WaterCourseElement * FindUpStreamWaterCourse(DistributedElement * const Dew, WaterCourseElement * const WaterWay, int * flowIndex, int nRows, int nCols, int noData, int row, int col);
void WriteWaterCourseFlow(DistributedElement * const Dew, WaterCourseElement * const WaterWay, WaterCourseElement ** Outlet, int landIndex, int flowIndex, int numOut, int nRows, int nCols, int noData, ofstream &fout);
void TraverseWaterCourse(WaterCourseElement * const thisWaterCourse);
void TraverseLandScape(DistributedElement * const thisElement, int waterCourseElementIdentification);
void WriteWaterCourseLandScapeElementsInformation(DistributedElement * const Dew, int landIndex, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);

int main(int argc, char *argv[])
{
  std::cout << "\n\n Preprocessing for Distributed Element Water Model \n\n";

  char fileName[80];
  char buffer[1024];
  char ch;
  int i,j,k;
  int numOut;
  int landIndex;
  int flowIndex;
  int nRows, nCols, noData;
  int numCatch, numCatchUp, numCatchOut;
  int stId1, stId2;
  double value;
  double xllCorner, yllCorner, cellSize;
  double correction, maxBas;

  if (argc != 2) {
    cout << " " << argv[0] << "  <control file name>\n\n";
    exit(1);
  }
  ifstream fileControl(argv[1]);
  if (fileControl == NULL) {
    cout << " Error opening file " << argv[1] << endl << endl;
    exit (1);
  }

  /*  cout << " Output file: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ofstream fout(fileName);  // Open for writing

  // Read landscape element file and generate landscape element objects  
  /*  cout << " File with geographical analysis area: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finGeo(fileName);  // Open for reading
  if (finGeo == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  } 
  finGeo >> buffer >> nCols;
  finGeo >> buffer >> nRows;
  finGeo >> buffer >> xllCorner;
  finGeo >> buffer >> yllCorner;
  finGeo >> buffer >> cellSize;
  finGeo >> buffer >> noData;
  //  cout << endl << nCols << endl;
  //  cout << nRows << endl;
  //  cout << xllCorner << endl;
  //  cout << yllCorner << endl;
  //  cout << cellSize << endl;
  //  cout << buffer << "    " << noData << endl;
  landIndex = 0;
  DistributedElement * Dew = new DistributedElement [nRows*nCols];
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      Dew[ELEMENT(i,j)].SetGeoIndex(ELEMENT(i,j));
      Dew[ELEMENT(i,j)].SetFlowIndex(noData);
      finGeo >> value;
      if (value!=noData) {
        Dew[ELEMENT(i,j)].SetLandIndex(landIndex++);
      } else {
        Dew[ELEMENT(i,j)].SetLandIndex(noData);
      }
    }
  }
  finGeo.close();
  // Read end

  // Read file with watershed outlets
  /*  cout << " File with watershed outlets: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finOut(fileName);  // Open for reading
  if (finOut == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  } 
  finOut.ignore(100,':');
  finOut >> numOut;
  int * outRow = new int[numOut];
  int * outCol = new int[numOut];
  for (i=0; i<numOut; i++) {
    finOut >> outRow[i] >> outCol[i];
    //    cout << outRow[i] << " " << outCol[i] << endl;
  }
  finOut.close();
  // Read end
  
  // Read watercourse elements identification, flow direction for landscape elements and flow direction for river elements
  ReadFlowDirection(Dew, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);
    
  // Find flow direction for landscape elements
  FindUpLandFlow(Dew, nRows, nCols, noData);
  
  // Generate watercourse element objects  
  WaterCourseElement * WaterWay = new WaterCourseElement [nRows*nCols];
  for (i=0; i<nRows*nCols; i++) {
    WaterWay[i].SetGeoIndex(Dew[i].GetGeoIndex());
    WaterWay[i].SetFlowIndex(noData);
  }
  
  // Watershed outlets and drainage network
  flowIndex = 0;
  WaterCourseElement ** Outlet = new WaterCourseElement * [numOut];
  for (i=0; i<numOut; i++) { 
    Outlet[i] = FindUpStreamWaterCourse(Dew, WaterWay, &flowIndex, nRows, nCols, noData, outRow[i], outCol[i]);
    if (!Outlet[i]) cout << "\n\n" << " Outlet[" << i << "] not found!" << "\n\n";
  }
  
  fout << endl;
  fout << "waterCourseElementIdentification orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetWaterCourseElementIdentification() << "  ";
    }
    fout << endl;
  }
  fout << endl;
  
  fout << "lakeNumber orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetLakeNumber() << "  ";
    }
    fout << endl;
  }
  fout << endl;
   
  fout << "flowLandValue orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetFlowLandValue() << "  ";
    }
    fout << endl;
  }
  fout << endl;
  
  fout << "flowRiverValue orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetFlowRiverValue() << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout <<"upLandFlow: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (Dew[ELEMENT(i,j)].GetNumUpLand()>0) {
        fout.width(10);
        fout << Dew[ELEMENT(i,j)].GetGeoIndex();
        for (k=0; k<Dew[ELEMENT(i,j)].GetNumUpLand(); k++) {
          fout.width(10);
          fout << Dew[ELEMENT(i,j)].GetUpLandFlow(k)->GetGeoIndex();
        }
        fout << endl;
      }
    }
  }
  fout << endl;
  
  fout <<"upStreamElement: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (Dew[ELEMENT(i,j)].GetFlowRiverValue()>=0) {
        fout.width(10);
        fout << WaterWay[ELEMENT(i,j)].GetGeoIndex();
        fout.width(10);
        fout << WaterWay[ELEMENT(i,j)].GetLandScapeElement()->GetGeoIndex();
        if (WaterWay[ELEMENT(i,j)].GetNumUpStream()>0) {
          for (k=0; k<WaterWay[ELEMENT(i,j)].GetNumUpStream(); k++) {
            fout << "   U";
            fout.width(10); 
            fout << WaterWay[ELEMENT(i,j)].GetUpStream(k)->GetGeoIndex();
            fout.width(10); 
            fout << WaterWay[ELEMENT(i,j)].GetUpStream(k)->GetFlowIndex();
          }
        }
        fout << endl;
      }
    }
  }
  fout << endl;
  cout << endl;

  // Write watercourse hierarchy elements identification
  WriteWaterCourseFlow(Dew, WaterWay, Outlet, landIndex, flowIndex, numOut, nRows, nCols, noData, fout);

  // Write watercourse/sub-catchment information for landscape elements 
  for (i=0; i<numOut; i++) { 
    TraverseWaterCourse(Outlet[i]);
  }
  WriteWaterCourseLandScapeElementsInformation(Dew, landIndex, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
  
  fout.close();
  fileControl.close();
  delete [] WaterWay;
  delete [] Outlet;
  delete [] outRow;
  delete [] outCol;
  delete [] Dew;
  return 0;
}


void GetRowCol(int elementNo, int nCols, int &row, int &col) 
{
  row = elementNo/nCols; 
  col = elementNo%nCols; 
}


int ControlFlowDirection(int flowDirection, int noData)
{
  return (flowDirection==1 || flowDirection==2 || flowDirection==4 || flowDirection==8 ||
          flowDirection==16 || flowDirection==32 || flowDirection==64 || flowDirection==128 ||
          flowDirection==noData);
}


void ReadFlowDirection(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[1024];
  int i,j,k;
  int nRo, nCo, noDa;
  int value;
  double xllC, yllC, cellS;

  // Read watercourse hierarchy elements identification
  /*  cout << endl << " File with watercourse hierarchy elements identification: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finWaterCourseElements(fileName);  // Open for reading
  if (finWaterCourseElements == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finWaterCourseElements >> buffer >> nCo;
  finWaterCourseElements >> buffer >> nRo;
  finWaterCourseElements >> buffer >> xllC;
  finWaterCourseElements >> buffer >> yllC;
  finWaterCourseElements >> buffer >> cellS;
  finWaterCourseElements >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for landscape elements!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finWaterCourseElements >> value;
    Dew[i].SetWaterCourseElementIdentification(value);
    if (Dew[i].GetWaterCourseElementIdentification()!=noData) {
      GetRowCol(i, nCols, j, k);
      //      fout << " WaterCourseElementIdentification = " <<  Dew[i].GetWaterCourseElementIdentification()<< "    Element : " << j << "," << k << "\n";
    }
  }
  finWaterCourseElements.close();
  // Read end

  // Read lake element numbers
  /*  cout << " File with lake element numbers: ";
      cin >> fileName;
      cout << endl;*/ 
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finLakeNumber(fileName);  // Open for reading
  if (finLakeNumber == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    //exit(1);
  }
  finLakeNumber >> buffer >> nCo;
  finLakeNumber >> buffer >> nRo;
  finLakeNumber >> buffer >> xllC;
  finLakeNumber >> buffer >> yllC;
  finLakeNumber >> buffer >> cellS;
  finLakeNumber >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for lake elements!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    //exit(1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finLakeNumber >> value;
    Dew[i].SetLakeNumber(value);
  }
  finLakeNumber.close();
  // Read end

  // Read flow direction grid for landscape elements
  /*  cout << endl << " File with flow direction grid for landscape elements: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finFlowLand(fileName);  // Open for reading
  if (finFlowLand == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finFlowLand >> buffer >> nCo;
  finFlowLand >> buffer >> nRo;
  finFlowLand >> buffer >> xllC;
  finFlowLand >> buffer >> yllC;
  finFlowLand >> buffer >> cellS;
  finFlowLand >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for landscape elements!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finFlowLand >> value;
    if (!ControlFlowDirection(value, noData)) {
      cout << "\n\n Error in flow direction for landscape elements!\n\n";
      cout << " value     " << value << endl;
      cout << " noData    " << noData << endl << endl;
      exit (1);
    }
    Dew[i].SetFlowLandValue(value);
    if (Dew[i].GetWaterCourseElementIdentification()!=noData) {
      Dew[i].SetFlowRiverValue(value);
      GetRowCol(i, nCols, j, k);
      //      fout << " FlowRiverValue = " <<  Dew[i].GetFlowRiverValue()<< "    Element : " << j << "," << k << "\n";
    }
    else Dew[i].SetFlowRiverValue(noData);
  }
  finFlowLand.close();
  // Read end

  // Read flow direction grid for river elements 
  /*
  cout << " File with flow direction grid for river elements: ";
  cin >> fileName;
  cout << endl;
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(1024,'\n');
  ifstream finFlowRiver(fileName);  // Open for reading
  if (finFlowRiver == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finFlowRiver >> buffer >> nCo;
  finFlowRiver >> buffer >> nRo;
  finFlowRiver >> buffer >> xllC;
  finFlowRiver >> buffer >> yllC;
  finFlowRiver >> buffer >> cellS;
  finFlowRiver >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for river elements!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finFlowRiver >> value;
    if (!ControlFlowDirection(value, noData)) {
      cout << "\n\n Error in flow direction for river elements!\n\n";
      cout << " value     " << value << endl;
      cout << " noData    " << noData << endl << endl;
      exit (1);
    }
    Dew[i].SetFlowRiverValue(value);
  }
  finFlowRiver.close();
  */
  // Read end

}


void FindUpLandFlow(DistributedElement * const Dew, int nRows, int nCols, int noData)
{
  int i,j;
  int nUp;

  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (Dew[ELEMENT(i,j)].GetLandIndex()!=noData) {
        nUp=0;
        if ((j<nCols-1) 
            && (Dew[ELEMENT(i,j+1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i,j+1)].GetFlowRiverValue()==noData)
            && (Dew[ELEMENT(i,j+1)].GetFlowLandValue()==16)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i,j+1)]);
          nUp++;
        }
        if ((i<nRows-1 && j<nCols-1) 
            && (Dew[ELEMENT(i+1,j+1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i+1,j+1)].GetFlowRiverValue()==noData)
            && (Dew[ELEMENT(i+1,j+1)].GetFlowLandValue()==32)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i+1,j+1)]);
          nUp++;
        }
        if ((i<nRows-1)
            && (Dew[ELEMENT(i+1,j)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i+1,j)].GetFlowRiverValue()==noData)
            && (Dew[ELEMENT(i+1,j)].GetFlowLandValue()==64)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i+1,j)]);
          nUp++;
        }
        if ((i<nRows-1) && (j>0) 
            && (Dew[ELEMENT(i+1,j-1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i+1,j-1)].GetFlowRiverValue()==noData)
            && (Dew[ELEMENT(i+1,j-1)].GetFlowLandValue()==128)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i+1,j-1)]);
          nUp++;
        }
        if ((j>0) 
            && (Dew[ELEMENT(i,j-1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i,j-1)].GetFlowRiverValue()==noData)
            && (Dew[ELEMENT(i,j-1)].GetFlowLandValue()==1)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i,j-1)]);
          nUp++;
        }
        if ((i>0 && j>0) 
            && (Dew[ELEMENT(i-1,j-1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i-1,j-1)].GetFlowRiverValue()==noData)
            && (Dew[ELEMENT(i-1,j-1)].GetFlowLandValue()==2)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i-1,j-1)]);
          nUp++;
        }
        if ((i>0)
            && (Dew[ELEMENT(i-1,j)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i-1,j)].GetFlowRiverValue()==noData)
            && (Dew[ELEMENT(i-1,j)].GetFlowLandValue()==4)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i-1,j)]);
          nUp++;
        }
        if ((i>0 && j<nCols-1)
            && (Dew[ELEMENT(i-1,j+1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i-1,j+1)].GetFlowRiverValue()==noData)
            && (Dew[ELEMENT(i-1,j+1)].GetFlowLandValue()==8)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i-1,j+1)]);
          nUp++;
        }
        if (nUp > 7) {
          cout << "\n\n Error in function FindUpLandFlow!\n\n";
          cout << "\n nUp: " << nUp << endl << endl;
          exit (1);
        }      
        Dew[ELEMENT(i,j)].SetNumUpLand(nUp);
      }
    }
  }
  return;
}


WaterCourseElement * FindUpStreamWaterCourse(DistributedElement * const Dew, WaterCourseElement * const WaterWay, int * flowIndex, int nRows, int nCols, int noData, int row, int col)
{
  int nUp=0;
  WaterCourseElement * result=0;

  //  cout << row << " " << col << " " << Dew[ELEMENT(row,col)].GetFlowRiverValue() << endl;
  // Not watercourse element  
  if (Dew[ELEMENT(row,col)].GetFlowRiverValue()==noData) {
    return result;
  } 

  // Watercourse element
  WaterWay[ELEMENT(row,col)].SetLandScapeElement(&Dew[ELEMENT(row,col)]);
  WaterWay[ELEMENT(row,col)].SetFlowIndex(*flowIndex);
  Dew[ELEMENT(row,col)].SetFlowIndex(*flowIndex);
  (*flowIndex)++;

  if (col<nCols-1 && WaterWay[ELEMENT(row,col+1)].GetFlowIndex()==noData &&
      Dew[ELEMENT(row,col+1)].GetFlowRiverValue()==16 && Dew[ELEMENT(row,col)].GetFlowRiverValue()!=1) {
    result = FindUpStreamWaterCourse(Dew, WaterWay, flowIndex, nRows, nCols, noData, row, col+1);
    if (result) {
      WaterWay[ELEMENT(row,col)].SetUpStreamElement(nUp, result);
      nUp++;
    }
  }
  if (row<nRows-1 && col<nCols-1 && WaterWay[ELEMENT(row+1,col+1)].GetFlowIndex()==noData &&
      Dew[ELEMENT(row+1,col+1)].GetFlowRiverValue()==32 && Dew[ELEMENT(row,col)].GetFlowRiverValue()!=2) {
    result = FindUpStreamWaterCourse(Dew, WaterWay, flowIndex, nRows, nCols, noData, row+1, col+1);
    if (result) {
      WaterWay[ELEMENT(row,col)].SetUpStreamElement(nUp, result);
      nUp++;
    }
  }
  if (row<nRows-1 && WaterWay[ELEMENT(row+1,col)].GetFlowIndex()==noData &&
      Dew[ELEMENT(row+1,col)].GetFlowRiverValue()==64 && Dew[ELEMENT(row,col)].GetFlowRiverValue()!=4) {
    result = FindUpStreamWaterCourse(Dew, WaterWay, flowIndex, nRows, nCols, noData, row+1, col);
    if (result) {
      WaterWay[ELEMENT(row,col)].SetUpStreamElement(nUp, result);
      nUp++;
    }
  }
  if (row<nRows-1 && col>0 && WaterWay[ELEMENT(row+1,col-1)].GetFlowIndex()==noData &&
      Dew[ELEMENT(row+1,col-1)].GetFlowRiverValue()==128 && Dew[ELEMENT(row,col)].GetFlowRiverValue()!=8) {
    result = FindUpStreamWaterCourse(Dew, WaterWay, flowIndex, nRows, nCols, noData, row+1, col-1);
    if (result) {
      WaterWay[ELEMENT(row,col)].SetUpStreamElement(nUp, result);
      nUp++;
    }
  }
  if (col>0 && WaterWay[ELEMENT(row,col-1)].GetFlowIndex()==noData &&
      Dew[ELEMENT(row,col-1)].GetFlowRiverValue()==1 && Dew[ELEMENT(row,col)].GetFlowRiverValue()!=16) {
    result = FindUpStreamWaterCourse(Dew, WaterWay, flowIndex, nRows, nCols, noData, row, col-1);
    if (result) {
      WaterWay[ELEMENT(row,col)].SetUpStreamElement(nUp, result);
      nUp++;
    }
  }
  if (row>0 && col>0 && WaterWay[ELEMENT(row-1,col-1)].GetFlowIndex()==noData &&
      Dew[ELEMENT(row-1,col-1)].GetFlowRiverValue()==2 && Dew[ELEMENT(row,col)].GetFlowRiverValue()!=32) {
    result = FindUpStreamWaterCourse(Dew, WaterWay, flowIndex, nRows, nCols, noData, row-1, col-1);
    if (result) {
      WaterWay[ELEMENT(row,col)].SetUpStreamElement(nUp, result);
      nUp++;
    }
  }
  if (row>0 && WaterWay[ELEMENT(row-1,col)].GetFlowIndex()==noData &&
      Dew[ELEMENT(row-1,col)].GetFlowRiverValue()==4 && Dew[ELEMENT(row,col)].GetFlowRiverValue()!=64) {
    result = FindUpStreamWaterCourse(Dew, WaterWay, flowIndex, nRows, nCols, noData, row-1, col);
    if (result) {
      WaterWay[ELEMENT(row,col)].SetUpStreamElement(nUp, result);
      nUp++;
    }
  }
  if (row>0 && col<nCols-1 && WaterWay[ELEMENT(row-1,col+1)].GetFlowIndex()==noData &&
      Dew[ELEMENT(row-1,col+1)].GetFlowRiverValue()==8 && Dew[ELEMENT(row,col)].GetFlowRiverValue()!=128) {
    result = FindUpStreamWaterCourse(Dew, WaterWay, flowIndex, nRows, nCols, noData, row-1, col+1);
    if (result) {
      WaterWay[ELEMENT(row,col)].SetUpStreamElement(nUp, result);
      nUp++;
    }
  }

  if (nUp > 7) {
    cout << "\n\n Error in function FindUpStreamWaterCourse!\n\n";
    cout << "\n nUp: " << nUp << endl << endl;
    exit (1);
  } else {
    WaterWay[ELEMENT(row,col)].SetNumUpStream(nUp);
  }
  return &WaterWay[ELEMENT(row,col)];
}


void WriteWaterCourseFlow(DistributedElement * const Dew, WaterCourseElement * const WaterWay, WaterCourseElement ** Outlet, int landIndex, int flowIndex, int numOut, int nRows, int nCols, int noData, ofstream &fout)
{
  int i, j, numLand, numWaterCourse;

  numWaterCourse=0;
  ofstream waterCourseHierarchyOut("dew_watercourse_hierarchy.txt");  // Open for writing
  waterCourseHierarchyOut << "# Number of watercourse hierarchy elements:  " << flowIndex << endl;
  for (i=0; i<nRows*nCols; i++) {
    if (WaterWay[i].GetFlowIndex()!=noData) {
      WaterWay[i].SetNumWaterCourse(numWaterCourse);
      //      waterCourseHierarchyOut << WaterWay[i].GetFlowIndex() << " : " << WaterWay[i].GetGeoIndex() << "  1.0 " << endl;
      waterCourseHierarchyOut << WaterWay[i].GetNumWaterCourse() << " : " << WaterWay[i].GetLandScapeElement()->GetWaterCourseElementIdentification() << "  " << WaterWay[i].GetLandScapeElement()->GetLakeNumber() << "  1.0  " << endl;
      numWaterCourse++;
    }
  }
  waterCourseHierarchyOut << "# Number of watercourse outlets:  " << numOut << endl;
  for (i=0; i<numOut; i++)
    waterCourseHierarchyOut << Outlet[i]->GetNumWaterCourse() << endl;
  waterCourseHierarchyOut << "# Hierarchy of watercourse elements " << endl;
  for (i=0; i<nRows*nCols; i++) {
    if (WaterWay[i].GetNumUpStream()>0) {
      waterCourseHierarchyOut << WaterWay[i].GetNumWaterCourse() << "  " << WaterWay[i].GetNumUpStream() << " :  ";
      for (j=0; j<WaterWay[i].GetNumUpStream(); j++)
        waterCourseHierarchyOut << WaterWay[i].GetUpStream(j)->GetNumWaterCourse() << "  ";
      waterCourseHierarchyOut << endl;
    }
  }
  waterCourseHierarchyOut.close();
}


void TraverseWaterCourse(WaterCourseElement * const thisWaterCourse)
{
  int i;
  for (i=0; i<thisWaterCourse->GetNumUpStream(); i++) { 
    TraverseWaterCourse(thisWaterCourse->GetUpStream(i));
  } 
  //  cout << " 1 " << thisWaterCourse->GetLandScapeElement()->GetWaterCourseElementIdentification() << endl;
  TraverseLandScape(thisWaterCourse->GetLandScapeElement(), thisWaterCourse->GetLandScapeElement()->GetWaterCourseElementIdentification());
}


void TraverseLandScape(DistributedElement * const thisElement, int waterCourseElementIdentification)
{
  int i;
  //  cout << " 2 " << waterCourseElementIdentification << endl;
  thisElement->SetWaterCourseElementIdentification(waterCourseElementIdentification);
  for (i=0; i<thisElement->GetNumUpLand(); i++) {
    TraverseLandScape(thisElement->GetUpLandFlow(i), waterCourseElementIdentification);
  }
 }


void WriteWaterCourseLandScapeElementsInformation(DistributedElement * const Dew, int landIndex, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
  int i, j;

  // Landscape elements  
  ofstream landScapeWaterCourseOut("dew_watercourse_landscape.txt");  // Open for writing
  landScapeWaterCourseOut << "ncols         " << nCols << endl;
  landScapeWaterCourseOut << "nrows         " << nRows << endl;
  landScapeWaterCourseOut << "xllcorner     " << xllCorner << endl;
  landScapeWaterCourseOut << "yllcorner     " << yllCorner << endl;
  landScapeWaterCourseOut << "cellsize      " << cellSize << endl;
  landScapeWaterCourseOut << "NODATA_value  " << noData << endl;
  //  for (i=0; i<nRows*nCols; i++) {
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (Dew[ELEMENT(i,j)].GetLandIndex()!=noData) {
        //      landScapeWaterCourseOut.width(10); 
        landScapeWaterCourseOut << Dew[ELEMENT(i,j)].GetWaterCourseElementIdentification() << endl;
      }
      else {
        //      landScapeWaterCourseOut.width(10); 
        landScapeWaterCourseOut << noData << endl;
      }
    }
    //    landScapeWaterCourseOut << endl ;
  }
  landScapeWaterCourseOut.close();
}
