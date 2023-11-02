/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments,          *
 *  river networks, lakes and landscape elements.                                       *
 *  Preprocessing.                                                                      *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "date_time.h"
#include "classPreDew.h"

void SetCommonParameters(ParametersCommon * const ParCommonStore, ifstream &fileControl, ofstream &fout);
void GetRowCol(int elementNo, int nCols, int &row, int &col); 
void ReadFlowDirection(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
void ReadLandUse(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
int ControlFlowDirection(int, int);
void FindUpLandFlow(DistributedElement * const Dew, int nRows, int nCols, int noData);
void WriteLandScapeElements(DistributedElement * const Dew, ParametersCommon * const ParCommonStore,
                            MeteorologicalStations * const MetStations, int landIndex, 
                            int preSetModelStructure, int nRows, int nCols, 
                            int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void FindCatchmentIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, int numWatc, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
void WriteCatchmentIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, int numWatc, ofstream &fout);
void FindWaterCourseIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, int numWatc, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout);
void WriteWaterCourseIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, int numWatc, ofstream &fout);
void WriteLandScapeAndWaterCourseFlow(DistributedElement * const Dew, int nRows, int nCols, int noData, ofstream &fout);
void DistanceSort(DistributedElement  * const distElement, ParametersCommon * const ParCommonStore, 
                  MeteorologicalStations * const MetStations,
                  int * precStations, int * tempStations, double * precWeights, double * tempWeights);

int main(int argc, char *argv[])
{
  std::cout << "\n\n Preprocessing for Distributed Element Water Model \n\n";

  char fileName[80];
  char buffer[256];
  char hierarchy = 'H';
  char ch;
  int i,j,k;
  int landIndex;
  int nRows, nCols, noData;
  int numWatc, numWatcUp, numWatcOut;
  int waterCourseId;
  int preSetModelStructure = -1;
  int lakeNumber;
  double value;
  double correction;
  double xllCorner, yllCorner, cellSize;

  DateTime nowDate;
  nowDate.now();
  DateTime finalDate(finalYear,finalMonth,finalDay,finalHour,finalMinute,0);
  //  cout<<" "<<finalDate.getYear()<<" "<<finalDate.getMonth()<<" "<<finalDate.getDay()<<" "<<finalDate.getHour()<<" "<<finalDate.getMinute()<<endl;
/*  if (nowDate > finalDate) {
    strcpy(buffer,"del ");
    strcat(buffer,argv[0]);
    strcat(buffer,"*\n");
    strcat(buffer,"del predew.bat\n");
    ofstream fexit("abc_987.bat");
    fexit << "@echo off\n" << buffer << endl;
    fexit.close();
    system("copy abc_987.bat predew.bat /Y >NUL");
    system("del abc_987.bat >NUL");
    cout << "\n *\n *  Program is terminated !  \n *\n\n\n";
    exit(1);
    }*/
  if (nowDate > finalDate) {
    remove(argv[0]);
    cout << "\n *\n *  Program is terminated !  \n *\n\n\n";
    exit(1);
  }

  if (argc != 2) {
    cout << " " << argv[0] << "  <control file name>\n\n";
    exit(1);
  }
  ifstream fileControl(argv[1]);
  if (fileControl == NULL) {
    cout << " Error opening file " << argv[1] << endl << endl;
    exit (1);
  }

  /*  while (preSetModelStructure < 0 || preSetModelStructure >= numberModelStructures) {
      cout << " Model structure, HBV (0) or KinematicWave (1): ";
      cin >> preSetModelStructure;
      cout << endl;
      }*/   
  fileControl.ignore(100,':');
  fileControl >> preSetModelStructure;
  fileControl.ignore(256,'\n');
  if (preSetModelStructure < 0 || preSetModelStructure >= numberModelStructures) {
    cout << "\n Model structure, HBV (0) or KinematicWave (1): \n\n";
    exit (1);
  }

  /*  while (hierarchy != 'N' && hierarchy != 'C' && hierarchy != 'n' && hierarchy != 'c') {
      cout << " Watercourse hierarchy description, river/lake network(N) or nested catchments(C): ";
      cin >> hierarchy;
      cout << endl;
      }*/
  fileControl.ignore(100,':');
  fileControl >> hierarchy;
  fileControl.ignore(256,'\n');
  if (hierarchy != 'N' && hierarchy != 'C' && hierarchy != 'n' && hierarchy != 'c') {
    cout << "\n Landscape elements hierarchy, flow direction network(N) or nested catchments(C) \n\n";
    exit (1);
  }

  /*  cout << " Output file: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ofstream fout(fileName);  // Open for writing

  // Object for storing meteorological station information
  MeteorologicalStations * MetStations = new MeteorologicalStations;
  MetStations->SetMeteorologicalStations(fileControl, fout);

  // Read common parameters file and set parameter values
  ParametersCommon * ParCommonStore = new ParametersCommon;
  SetCommonParameters(ParCommonStore, fileControl, fout);
  // End read common parameters

  // Read landscape element file and generate landscape element objects  
  /*  cout << " File with geographical analysis area: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
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
  cout << endl << nCols << endl;
  cout << nRows << endl;
  cout << xllCorner << endl;
  cout << yllCorner << endl;
  cout << cellSize << endl;
  cout << buffer << "    " << noData << endl;
  landIndex = 0;
  DistributedElement * Dew = new DistributedElement [nRows*nCols];
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      Dew[ELEMENT(i,j)].SetGeoIndex(ELEMENT(i,j));
      finGeo >> value;
      if (value!=noData) {
        Dew[ELEMENT(i,j)].SetLandIndex(landIndex++);
      } else {
        Dew[ELEMENT(i,j)].SetLandIndex(noData);
      }
      Dew[ELEMENT(i,j)].SetXCoord(xllCorner + (j+0.5)*cellSize);
      Dew[ELEMENT(i,j)].SetYCoord(yllCorner + (nRows-i-0.5)*cellSize);
    }
  }
  finGeo.close();
  // Read end

  // Read land use for landscape elements
  ReadLandUse(Dew, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);

  // Write landscape element information to output file
  /*  fout << endl << "geoIndex orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      if (Dew[ELEMENT(i,j)].GetLandIndex() != noData)
        fout << Dew[ELEMENT(i,j)].GetGeoIndex() << "  ";
      else
        fout << noData << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout << "landIndex orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetLandIndex() << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout << "rader og kolonner orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout << "   " << i << "," << j;
    }
    fout << endl;
  }
  fout << endl;

  fout << "elevation orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetElevation() << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout << "slopeLength orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetSlopeLength() << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout << "slopeAngle orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetSlopeAngle() << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout << "lakePercent orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetLakePercent() << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout << "forestPercent orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetForestPercent() << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout << "bogPercent orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetBogPercent() << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout << "glacierPercent orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetGlacierPercent() << "  ";
    }
    fout << endl;
  }
  fout << endl;

  fout << "openLandPercent orientert nord-sør: \n";
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      fout.width(5);
      fout << Dew[ELEMENT(i,j)].GetOpenLandPercent() << "  ";
    }
    fout << endl;
  }
  fout << endl;*/

  // Open file with watercourse/catchment hierarchy
  /*  cout << " File with watercourse/catchment hierarchy: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream fileWCo(fileName);
  if (fileWCo == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  
  // Read watercourse/catchment information and generate watercourse objects  
  fileWCo.ignore(100,':');
  fileWCo >> numWatc;
  cout << "\n # Number of watercourses/catchments " << numWatc << endl;
  WaterCourse *WaterElement = new WaterCourse [numWatc];
  for (i=0; i<numWatc; i++) {
    fileWCo >> j >> ch >> waterCourseId >> lakeNumber >> correction;
    if (j != i) {
      cout << endl << "Error reading file " << fileName << "\t" << i << "\t" << j << endl;
      exit (1);
    }
    fileWCo.ignore(256,'\n');
    WaterElement[i].SetWaterCourseIndex(i);
    WaterElement[i].SetIdentifier(waterCourseId);
    WaterElement[i].SetCorrection(correction);
    cout << " Watercourse/catchment index " << i << "  " << "Watercourse/catchment identifier " << WaterElement[i].GetIdentifier() << endl;
  }

  // Watercourse outlets
  fileWCo.ignore(100,':');
  fileWCo >> numWatcOut;
  cout << "\n # Number of watercourse/catchment outlets " << numWatcOut << endl;
  WaterCourse ** Outlet = new WaterCourse * [numWatcOut];
  for (i=0; i<numWatcOut; i++) {
    fileWCo >> j;
    Outlet[i] = &WaterElement[j];
    cout << " Outlet no. " << i << "\t" << " Watercourse/catchment no. " << j << "\t" << endl;
    fileWCo.ignore(256,'\n');
  }

  // Hierarchy of watercourses
  fileWCo.getline(buffer, 256);
  cout << "\n " << buffer << endl;
  while (fileWCo >> i) {
    fileWCo >> numWatcUp;
    WaterElement[i].SetNumUpStream(numWatcUp);
    fileWCo.ignore(100,':');
    cout << " Downstream, watercourse/catchment no.  " << i << "    Identifier  " << WaterElement[i].GetIdentifier() << endl; 
    cout << " No. of upstream watercourses " << numWatcUp << endl;
    k = 0;
    while (fileWCo.peek() != '\n') {
      fileWCo >> j;
      cout << "\t" << "Upstream, watercourse/catchment no. " << j ;
      WaterElement[i].SetUpStream(k, &WaterElement[j]);
      cout  << "\t" << "UpStream[" << k << "]" << "    Identifier  " << WaterElement[i].GetUpStream(k)->GetIdentifier() << endl;
      while (fileWCo.peek() == ' ') fileWCo.ignore(1,' ');
      k++;
    }
    fileWCo.ignore(256,'\n');
    if (numWatcUp!=k) {
      cout << endl << " Error in number of upstream pointers for watercourse/catchment no. " << i << endl << endl;
      exit (1);
    } 
  }
  fileWCo.close();
  // Read end

  // Read flow direction for landscape elements
  ReadFlowDirection(Dew, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);
    
  // Write landscape element information to input file to dew
  WriteLandScapeElements(Dew, ParCommonStore, MetStations, landIndex, preSetModelStructure, 
                         nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
  
  if (hierarchy == 'N' || hierarchy == 'n') {
    // Find watercourse identifiers for landscape elements and 
    // connect landscape elements to watercourse elements
    FindWaterCourseIdentifier(Dew, WaterElement, numWatc, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout); 
    
    // Write information about watercourse elements and landscape elements to input file to dew.cpp
    WriteWaterCourseIdentifier(Dew, WaterElement, numWatc, fout);
    
    // Find flow direction for landscape elements
    FindUpLandFlow(Dew, nRows, nCols, noData);
    
    // Write flow direction for landscape elements
    WriteLandScapeAndWaterCourseFlow(Dew, nRows, nCols, noData, fout);
  } 
  
  else if (hierarchy == 'C' || hierarchy == 'c') {
    // Find catchment identifiers for landscape elements and 
    // connect landscape elements to catchment outlets
    FindCatchmentIdentifier(Dew, WaterElement, numWatc, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fileControl, fout);
    
    // Write information about catchment outlets and landscape elements to input file to dew.cpp
    WriteCatchmentIdentifier(Dew, WaterElement, numWatc, fout);
  }
  
  delete [] Outlet;
  delete [] WaterElement;
  
  fout.close();
  fileControl.close();
  delete [] Dew;
  return 0;
}


void SetCommonParameters(ParametersCommon * const ParCommonStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  int numPrec, numTemp;

  /*  cout << " File with common parameters: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finCommonPar(fileName);  // Open for reading
  if (finCommonPar == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finCommonPar.ignore(256,'\n');
  finCommonPar.ignore(100,':'); finCommonPar >> numPrec; 
  finCommonPar.ignore(100,':'); finCommonPar >> numTemp; 
  ParCommonStore->SetNUM_PREC_SERIES(numPrec);
  ParCommonStore->SetNUM_TEMP_SERIES(numTemp);
  finCommonPar.close();
  fout << "Common parameters: \n";
  fout << ParCommonStore->GetNUM_PREC_SERIES() << endl;
  fout << ParCommonStore->GetNUM_TEMP_SERIES() << endl;
  fout << endl << endl;
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


void ReadLandUse(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i,j,k;
  int nRo, nCo, noDa;
  double value, notAssigned;
  double xllC, yllC, cellS;

  /*  cout << endl;
  for (i=0; i<nRows*nCols; i++) {
    cout << Dew[i].GetGeoIndex() << "  ";
  }
  cout << endl << endl;

  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      cout << Dew[ELEMENT(i,j)].GetGeoIndex() << "  ";
    }
    cout << endl;
  }
  cout << endl;

  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      cout << ELEMENT(i,j) << "  ";
    }
  }
  cout << endl << endl;

  for (k=0; k<nRows*nCols; k++) {
    GetRowCol(k, nCols, i, j);
    if (k>0 && j%nCols==0) cout << endl;
    cout << " " << i << "," << j;
  }
  cout << endl << endl; */

  // Area of grid cells 
  for (i=0; i<nRows*nCols; i++) {
    Dew[i].SetArea(cellSize*cellSize);
  }

  // Read elevation of grid cells
  /*  cout << " File with grid cell elevations: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finElevation(fileName);  // Open for reading
  if (finElevation == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finElevation >> buffer >> nCo;
  finElevation >> buffer >> nRo;
  finElevation >> buffer >> xllC;
  finElevation >> buffer >> yllC;
  finElevation >> buffer >> cellS;
  finElevation >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for elevation!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finElevation >> value;
    if (value<0.0) value=0.0;
    Dew[i].SetElevation(value);
  }
  finElevation.close();
  // Read end

  // Read slope length of grid cells
  /*  cout << " File with slope lengths: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finSlopeLength(fileName);  // Open for reading
  if (finSlopeLength == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finSlopeLength >> buffer >> nCo;
  finSlopeLength >> buffer >> nRo;
  finSlopeLength >> buffer >> xllC;
  finSlopeLength >> buffer >> yllC;
  finSlopeLength >> buffer >> cellS;
  finSlopeLength >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for slope!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finSlopeLength >> value;
    if (value<0) value=0;
    Dew[i].SetSlopeLength(value);
  }
  finSlopeLength.close();
  // Read end

  // Read slope angle of grid cells
  /*  cout << " File with slope angles: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finSlopeAngle(fileName);  // Open for reading
  if (finSlopeAngle == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finSlopeAngle >> buffer >> nCo;
  finSlopeAngle >> buffer >> nRo;
  finSlopeAngle >> buffer >> xllC;
  finSlopeAngle >> buffer >> yllC;
  finSlopeAngle >> buffer >> cellS;
  finSlopeAngle >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for slope!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finSlopeAngle >> value;
    if (value<0) value=0;
    Dew[i].SetSlopeAngle(value);
  }
  finSlopeAngle.close();
  // Read end

  // Read slope aspect of grid cells
  /*  cout << " File with slope aspects: ";
    cin >> fileName;
    cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finAspect(fileName);  // Open for reading
  if (finAspect == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finAspect >> buffer >> nCo;
  finAspect >> buffer >> nRo;
  finAspect >> buffer >> xllC;
  finAspect >> buffer >> yllC;
  finAspect >> buffer >> cellS;
  finAspect >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for aspect!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finAspect >> value;
    if (value<0) value=0;
    Dew[i].SetAspect(value);
  }
  finAspect.close();
  // Read end

  // Read percentage of grid cells covered by lakes
  /*  cout << " File with lake percentage: ";
      cin >> fileName;
      cout << endl; */
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finLakePercent(fileName);  // Open for reading
  if (finLakePercent == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finLakePercent >> buffer >> nCo;
  finLakePercent >> buffer >> nRo;
  finLakePercent >> buffer >> xllC;
  finLakePercent >> buffer >> yllC;
  finLakePercent >> buffer >> cellS;
  finLakePercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for lake percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finLakePercent >> value;
    if (value<0.0) value=0.0;
    Dew[i].SetLakePercent(value);
  }
  finLakePercent.close();
  // Read end

  // Read percentage of grid cells covered by forest
  /*  cout << " File with forest percentage: ";
      cin >> fileName;
      cout << endl; */
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finForestPercent(fileName);  // Open for reading
  if (finForestPercent == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finForestPercent >> buffer >> nCo;
  finForestPercent >> buffer >> nRo;
  finForestPercent >> buffer >> xllC;
  finForestPercent >> buffer >> yllC;
  finForestPercent >> buffer >> cellS;
  finForestPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for forest percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finForestPercent >> value;
    if (value<0.0) value=0.0;
    Dew[i].SetForestPercent(value);
  }
  finForestPercent.close();
  // Read end

  // Read percentage of grid cells covered by bogs
  /*  cout << " File with bog percentage: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finBogPercent(fileName);  // Open for reading
  if (finBogPercent == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finBogPercent >> buffer >> nCo;
  finBogPercent >> buffer >> nRo;
  finBogPercent >> buffer >> xllC;
  finBogPercent >> buffer >> yllC;
  finBogPercent >> buffer >> cellS;
  finBogPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for bog percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finBogPercent >> value;
    if (value<0.0) value=0.0;
    Dew[i].SetBogPercent(value);
  }
  finBogPercent.close();
  // Read end

  // Read percentage of grid cells covered by glaciers
  /*  cout << " File with glacier percentage: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finGlacierPercent(fileName);  // Open for reading
  if (finGlacierPercent == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finGlacierPercent >> buffer >> nCo;
  finGlacierPercent >> buffer >> nRo;
  finGlacierPercent >> buffer >> xllC;
  finGlacierPercent >> buffer >> yllC;
  finGlacierPercent >> buffer >> cellS;
  finGlacierPercent >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finGlacierPercent >> value;
    if (value<0.0) value=0.0;
    Dew[i].SetGlacierPercent(value);
  }
  finGlacierPercent.close();
  // Read end

  // Read tree level for grid cells
  /*  cout << " File with tree levels: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finTreeLevel(fileName);  // Open for reading
  if (finTreeLevel == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finTreeLevel >> buffer >> nCo;
  finTreeLevel >> buffer >> nRo;
  finTreeLevel >> buffer >> xllC;
  finTreeLevel >> buffer >> yllC;
  finTreeLevel >> buffer >> cellS;
  finTreeLevel >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for glacier percent!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finTreeLevel >> value;
    if (value<0.0) value=0.0;
    Dew[i].SetTreeLevel(value);
  }
  finTreeLevel.close();
  // Read end

  // Assign land use classes based on potential tree level
  for (i=0; i<nRows*nCols; i++) {
    if (Dew[i].GetLandIndex()!=noData) { 
      notAssigned = 100.0 - Dew[i].GetLakePercent() - Dew[i].GetGlacierPercent() -
        Dew[i].GetForestPercent() - Dew[i].GetBogPercent();

      if (Dew[i].GetElevation() > Dew[i].GetTreeLevel()) {
        if (Dew[i].GetForestPercent() > 0.0) {
          if (notAssigned > 0.0) {
            Dew[i].SetAlpinePercent(Dew[i].GetForestPercent() + notAssigned*0.2);
            Dew[i].SetHeatherPercent(notAssigned*0.8);
            Dew[i].SetForestPercent(0.0);
          } 
          else {
            Dew[i].SetAlpinePercent(Dew[i].GetForestPercent());
            Dew[i].SetForestPercent(0.0);
          }
        }
        else if (Dew[i].GetElevation() <= Dew[i].GetTreeLevel()+100.0) {
          Dew[i].SetAlpinePercent(notAssigned*0.2);
          Dew[i].SetHeatherPercent(notAssigned*0.8);
        }
        else if (Dew[i].GetElevation() <= Dew[i].GetTreeLevel()+200.0) {
          Dew[i].SetHeatherPercent(notAssigned*0.5);
          Dew[i].SetBedrockPercent(notAssigned*0.5);
        }
        else { 
          Dew[i].SetBedrockPercent(notAssigned);
        }
      }
      else if (Dew[i].GetElevation() > Dew[i].GetTreeLevel()*0.9) {
        if (notAssigned > 0.0) {
          Dew[i].SetAlpinePercent(Dew[i].GetForestPercent() + notAssigned*0.2);
          Dew[i].SetHeatherPercent(notAssigned*0.8);
          Dew[i].SetForestPercent(0.0);
        }
        else {
          Dew[i].SetAlpinePercent(Dew[i].GetForestPercent()*0.5);
          Dew[i].SetForestPercent(Dew[i].GetForestPercent()*0.5);
        }
      }
      else if (Dew[i].GetElevation() > Dew[i].GetTreeLevel()*0.8) {
        if (notAssigned > 0.0) {
          Dew[i].SetOpenLandPercent(notAssigned);
        }
        else {
          Dew[i].SetAlpinePercent(Dew[i].GetForestPercent()*0.2);
          Dew[i].SetForestPercent(Dew[i].GetForestPercent()*0.8);
        }
      }
      else {
        if (notAssigned > 0.0) {
          Dew[i].SetOpenLandPercent(notAssigned);
        }
      }
    }
  }

  // Control and correct area fractions of land use classes
  for (i=0; i<nRows*nCols; i++) {
    if (Dew[i].GetLandIndex()!=noData) { 
      value=Dew[i].GetForestPercent()+Dew[i].GetAlpinePercent()+
        Dew[i].GetHeatherPercent()+Dew[i].GetBedrockPercent()+
        Dew[i].GetLakePercent()+Dew[i].GetGlacierPercent()+Dew[i].GetBogPercent();

      if (value > 100.0) {
        Dew[i].SetForestPercent(Dew[i].GetForestPercent()*100.0/value);
        Dew[i].SetAlpinePercent(Dew[i].GetAlpinePercent()*100.0/value);
        Dew[i].SetHeatherPercent(Dew[i].GetHeatherPercent()*100.0/value);
        Dew[i].SetBedrockPercent(Dew[i].GetBedrockPercent()*100.0/value);
        Dew[i].SetLakePercent(Dew[i].GetLakePercent()*100.0/value);
        Dew[i].SetGlacierPercent(Dew[i].GetGlacierPercent()*100.0/value);
        Dew[i].SetBogPercent(Dew[i].GetBogPercent()*100.0/value);

        value=Dew[i].GetForestPercent()+Dew[i].GetAlpinePercent()+
        Dew[i].GetHeatherPercent()+Dew[i].GetBedrockPercent()+
        Dew[i].GetLakePercent()+Dew[i].GetGlacierPercent()+Dew[i].GetBogPercent();
      }

      if (value > 100.0) {
        GetRowCol(i, nCols, j, k);
        cout << "\n Sum of land use percentages = " << value << "    Element : " << j << "," << k << "\n";
        Dew[i].SetOpenLandPercent(0.0);
      } else {
        Dew[i].SetOpenLandPercent(100.0-value);
      }
    }
  }
}


void WriteLandScapeElements(DistributedElement * const Dew, ParametersCommon * const ParCommonStore,
                            MeteorologicalStations * const MetStations, int landIndex, 
                            int preSetModelStructure, int nRows, int nCols, int noData, 
                            double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
  int i, k, n, maxIndex;
  double lakePer, glacPer, areaCorr[2], totalArea, tempArea, areaFraction[numberLandSurfaceClasses];
  MODEL_STRUCTURE modelStructure[numberModelStructures];
  LANDSURFACE tempLandSurf, landSurfType[numberLandSurfaceClasses];
  SOIL soilType[numberSoilClasses];
  int * precipitationStation = new int [ParCommonStore->GetNUM_PREC_SERIES()];
  int * temperatureStation = new int [ParCommonStore->GetNUM_TEMP_SERIES()];
  double * precipitationWeight = new double [ParCommonStore->GetNUM_PREC_SERIES()];
  double * temperatureWeight = new double [ParCommonStore->GetNUM_TEMP_SERIES()];

  for (k=0; k<numberModelStructures; k++) modelStructure[k]=MODEL_STRUCTURE(k);
  // Landscape elements  
  ofstream landScapeOut("dew_landscape.txt");  // Open for writing
  ofstream geoIndexOut("dew_grid_index.txt");  // Open for writing
  landScapeOut << "ncols         " << nCols << endl;
  landScapeOut << "nrows         " << nRows << endl;
  landScapeOut << "xllcorner     " << xllCorner << endl;
  landScapeOut << "yllcorner     " << yllCorner << endl;
  landScapeOut << "cellsize      " << cellSize << endl;
  landScapeOut << "NODATA_value  " << noData << endl;
  landScapeOut << "# Number of landscape elements :  " << landIndex << endl;
  geoIndexOut << " Index of grid cells relative to upper left corner of rectangle with no. of rows = " << nRows << " and no. of columns = " << nCols << endl;
  for (i=0; i<nRows*nCols; i++) {
    if (Dew[i].GetLandIndex()!=noData) {
      // Sort land surface types based on area
      for (k=0; k<numberLandSurfaceClasses; k++) landSurfType[k]=LANDSURFACE(k);
      for (k=0; k<numberLandSurfaceClasses; k++) areaFraction[k] = 0.0;
      areaFraction[0]=Dew[i].GetOpenLandPercent();      // Open land (meadows, agriculture)
      areaFraction[1]=Dew[i].GetBogPercent();           // Bogs
      areaFraction[2]=Dew[i].GetForestPercent();        // Forest
      areaFraction[3]=Dew[i].GetAlpinePercent();        // Alpine forest
      areaFraction[4]=Dew[i].GetHeatherPercent();       // Low mountain with vegetation
      areaFraction[5]=Dew[i].GetBedrockPercent();       // Exposed bedrock (high mountain)
      /*      for (k=0; k<numberLandSurfaceClasses; k++) {
              maxIndex=k;
              for (n=k+1; n<numberLandSurfaceClasses; n++)
              if (areaFraction[n]>areaFraction[maxIndex]) maxIndex=n;
              tempLandSurf=landSurfType[maxIndex];
              landSurfType[maxIndex]=landSurfType[k];
              landSurfType[k]=tempLandSurf;
              tempArea=areaFraction[maxIndex];
              areaFraction[maxIndex]=areaFraction[k];
              areaFraction[k]=tempArea; 
              }*/
      for (k=0; k<numberLandSurfaceClasses; k++) {
        maxIndex=k;
        for (n=k+1; n<numberLandSurfaceClasses; n++) {
          if (areaFraction[n] > areaFraction[maxIndex]) {
            tempLandSurf=landSurfType[maxIndex];
            landSurfType[maxIndex]=landSurfType[n];
            landSurfType[n]=tempLandSurf;
            tempArea=areaFraction[maxIndex];
            areaFraction[maxIndex]=areaFraction[n];
            areaFraction[n]=tempArea; 
          }
        }
      }
      for (k=0; k<numberLandSurfaceClasses; k++) soilType[k]=SOIL(landSurfType[k]);
      lakePer=Dew[i].GetLakePercent();
      glacPer=Dew[i].GetGlacierPercent();
      areaCorr[0]=areaFraction[0];
      areaCorr[1]=areaFraction[1];
      //      cout << lakePer << "    " << glacPer << "    "  << areaCorr[0] << "    "
      //           << areaCorr[1] << endl;
      totalArea = lakePer + glacPer + areaCorr[0] + areaCorr[1];
      if (totalArea!=100.0) {
        //      if (areaFraction[0]==0.0 && areaFraction[1]==0.0) {
        //        lakePer = lakePer*100/totalArea;
        //        glacPer = glacPer*100/totalArea;
        //        } else {
        areaCorr[0]=areaFraction[0]+areaFraction[0]*(100.0-totalArea)/(areaFraction[0]+areaFraction[1]);
        areaCorr[1]=areaFraction[1]+areaFraction[1]*(100.0-totalArea)/(areaFraction[0]+areaFraction[1]);
        //       }
        totalArea = lakePer + glacPer + areaCorr[0] + areaCorr[1];
        if (totalArea!=100.0) cout << endl << "totalArea    " << totalArea << endl;
      }

      /* Find meteorological stations */
      DistanceSort(&Dew[i], ParCommonStore, MetStations, 
                   precipitationStation, temperatureStation, precipitationWeight, temperatureWeight);

      landScapeOut.precision(4); landScapeOut.setf(ios::fixed); landScapeOut.setf(ios::showpoint);
      landScapeOut.width(10); landScapeOut << Dew[i].GetLandIndex() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetGeoIndex() << "  ";
      landScapeOut.width(10); landScapeOut << modelStructure[preSetModelStructure] << "  ";
      landScapeOut.width(15); landScapeOut << Dew[i].GetArea() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetElevation() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetSlopeLength() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetSlopeAngle() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetAspect() << "  ";
      landScapeOut.width(10); landScapeOut << Dew[i].GetFlowLandValue() << "  ";
      landScapeOut.width(10); landScapeOut << lakePer << "  ";
      landScapeOut.width(10); landScapeOut << glacPer << "  ";
      for (k=0; k<maximumNumberLandClasses; k++) {
        landScapeOut.width(10); landScapeOut << landSurfType[k] << "  ";
        landScapeOut.width(10); landScapeOut << soilType[k] << "  ";
        landScapeOut.width(10); landScapeOut << areaCorr[k] << "  ";
      }
      for (k=0; k<ParCommonStore->GetNUM_PREC_SERIES(); k++) {
        landScapeOut.width(10); landScapeOut << precipitationStation[k] << "  ";
        landScapeOut.width(10); landScapeOut << precipitationWeight[k] << "  ";
      }
      for (k=0; k<ParCommonStore->GetNUM_TEMP_SERIES(); k++) {
        landScapeOut.width(10); landScapeOut << temperatureStation[k] << "  ";
        landScapeOut.width(10); landScapeOut << temperatureWeight[k] << "  ";
      }
      landScapeOut<< endl;
      geoIndexOut.width(10); geoIndexOut << "  " << Dew[i].GetGeoIndex() << "  " << endl;
    }
  }
  landScapeOut << endl;
  landScapeOut.close();
  geoIndexOut.close();

  delete [] precipitationStation;
  delete [] temperatureStation;
  delete [] precipitationWeight;
  delete [] temperatureWeight;
}


void FindCatchmentIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, int numWatc, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i,j;
  int nRo, nCo, noDa;
  int waterCourseId;
  double xllC, yllC, cellS;
  bool waterCourseFound;
  bool noElement=false;
  DistributedElement * lastElement;
  
  // Read watercourse/catchment identifiers for landscape elements
  /*  cout << "\n File with watercourse/catchment identifiers: ";
      cin >> file;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finId(fileName);  // Open for reading
  if (finId == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finId >> buffer >> nCo;
  finId >> buffer >> nRo;
  finId >> buffer >> xllC;
  finId >> buffer >> yllC;
  finId >> buffer >> cellS;
  finId >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for watercourse/catchment identifier one!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  // Read end

  // Connect landscape elements to watercourse/catchment outlets
  for (i=0; i<nRows*nCols; i++) {
    finId >> waterCourseId;
    Dew[i].SetWaterCourseValue(waterCourseId);
    j=0;
    waterCourseFound = false;
    while (j<numWatc) {
      if (waterCourseId==WaterElement[j].GetIdentifier()) {
        if (!waterCourseFound) {
          waterCourseFound=true;
          if (WaterElement[j].GetLandScapeElement()) {
            lastElement = WaterElement[j].GetLandScapeElement();
            while (lastElement->GetNextElement()) lastElement = lastElement->GetNextElement();
            lastElement->SetNextElement(&Dew[i]);
          } else {
            WaterElement[j].SetLandScapeElement(&Dew[i]);
          }
          WaterElement[j].SetNumLandScape(WaterElement[j].GetNumLandScape()+1);
        } else {
          cout << "\n\n Error in watercourse/catchment identifiers!   ";
          cout << WaterElement[j].GetIdentifier() << endl << endl;
          exit (1);
        }
      } 
      j++;
    }
  }

  for (j=0; j<numWatc; j++) {
    if (!(WaterElement[j].GetLandScapeElement())) {
      cout << "\n No landscape elements for watercourse/catchment index " << j << "  "  << "Watercourse/catchment identifier " << WaterElement[j].GetIdentifier() << endl;
      noElement=true;
    }
  }
  if (noElement) {
    cout << "\n Program is terminated! " << endl;
    exit(1);
  }
    
  finId.close();
}


void WriteCatchmentIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, int numWatc, ofstream &fout)
{
  int i;
  DistributedElement * thisElement;

  // Catchment outlets
  ofstream waterCourseOut("dew_waterland.txt");  // Open for writing
  for (i=0; i<numWatc; i++) {
    if (WaterElement[i].GetLandScapeElement()) {
      waterCourseOut << "#  " << WaterElement[i].GetIdentifier() << "  #  " << WaterElement[i].GetNumLandScape() << endl;
      thisElement = WaterElement[i].GetLandScapeElement();
      while (thisElement) {
        waterCourseOut << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
        thisElement = thisElement->GetNextElement();
      }
    }
  }
  waterCourseOut << endl;
  waterCourseOut.close();
}


void ReadFlowDirection(DistributedElement * const Dew, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i;
  int nRo, nCo, noDa;
  int value;
  double xllC, yllC, cellS;

  // Read flow direction grid for landscape elements
  /*  cout << endl << "\n File with flow direction grid for landscape elements: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
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
  }
  finFlowLand.close();
  // Read end
}


void FindWaterCourseIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, int numWatc, int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i,j;
  int nRo, nCo, noDa;
  int waterCourseId;
  double xllC, yllC, cellS;
  bool waterCourseFound;
  bool noElement=false;
  DistributedElement * lastElement;

  // Read watercourse identifiers for landscape elements
  /*  cout << "\n File with watercourse identifiers: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finId(fileName);  // Open for reading
  if (finId == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  }
  finId >> buffer >> nCo;
  finId >> buffer >> nRo;
  finId >> buffer >> xllC;
  finId >> buffer >> yllC;
  finId >> buffer >> cellS;
  finId >> buffer >> noDa;
  if (nCo!=nCols || nRo!=nRows || xllC!=xllCorner || yllC!=yllCorner || cellS!=cellSize || noDa!=noData) {
    cout << "\n\n Error in grid header information for watercourse elements!\n\n";
    cout << " nCols        " << nCo << "\t  " << nCols << endl;
    cout << " nRows        " << nRo << "\t  " << nRows << endl;
    cout << " xllCorner    " << xllC << "\t  " << xllCorner << endl;
    cout << " yllCorner    " << yllC << "\t  " << yllCorner << endl;
    cout << " cellSize     " << cellS << "\t  " << cellSize << endl;
    cout << " noData       " << noDa << "\t  " << noData << endl << endl;
    exit (1);
  }
  for (i=0; i<nRows*nCols; i++) {
    finId >> waterCourseId;
    Dew[i].SetWaterCourseValue(waterCourseId);
    j=0;
    waterCourseFound = false;
    while (j<numWatc) {
      if (waterCourseId==WaterElement[j].GetIdentifier()) {
        if (!waterCourseFound) {
          waterCourseFound=true;
          if (WaterElement[j].GetLandScapeElement()) {
            lastElement = WaterElement[j].GetLandScapeElement();
            while (lastElement->GetNextElement()) lastElement = lastElement->GetNextElement();
            lastElement->SetNextElement(&Dew[i]);
          } else {
            WaterElement[j].SetLandScapeElement(&Dew[i]);
          }
          WaterElement[j].SetNumLandScape(WaterElement[j].GetNumLandScape()+1);
        } else {
          cout << "\n\n Error in watercourse/catchment identifiers!   ";
          cout << WaterElement[j].GetIdentifier() << endl << endl;
          exit (1);
        }
      } 
      j++;
    }
  }

  for (j=0; j<numWatc; j++) {
    if (!(WaterElement[j].GetLandScapeElement())) {
      cout << "\n No landscape elements for watercourse/catchment index " << j << "  "  << "Watercourse/catchment identifier " << WaterElement[j].GetIdentifier() << endl;
      noElement=true;
    }
  }
  if (noElement) {
    cout << "\n Program is terminated! " << endl;
    exit(1);
  }
  finId.close();
    
}


void WriteWaterCourseIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, int numWatc, ofstream &fout)
{
  int i;
  DistributedElement * thisElement;

  // Watercourse elements
  ofstream waterCourseOut("dew_waterland.txt");  // Open for writing
  for (i=0; i<numWatc; i++) {
    if (WaterElement[i].GetLandScapeElement()) {
      waterCourseOut << "#  " << WaterElement[i].GetIdentifier() << "  #  " << WaterElement[i].GetNumLandScape() << endl;
      thisElement = WaterElement[i].GetLandScapeElement();
      while (thisElement) {
        waterCourseOut << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
        thisElement = thisElement->GetNextElement();
      }
    }
  }
  waterCourseOut << endl;
  waterCourseOut.close();
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
            && (Dew[ELEMENT(i,j+1)].GetWaterCourseValue()==noData)
            && (Dew[ELEMENT(i,j+1)].GetFlowLandValue()==16)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i,j+1)]);
          nUp++;
        }
        if ((i<nRows-1 && j<nCols-1) 
            && (Dew[ELEMENT(i+1,j+1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i+1,j+1)].GetWaterCourseValue()==noData)
            && (Dew[ELEMENT(i+1,j+1)].GetFlowLandValue()==32)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i+1,j+1)]);
          nUp++;
        }
        if ((i<nRows-1)
            && (Dew[ELEMENT(i+1,j)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i+1,j)].GetWaterCourseValue()==noData)
            && (Dew[ELEMENT(i+1,j)].GetFlowLandValue()==64)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i+1,j)]);
          nUp++;
        }
        if ((i<nRows-1) && (j>0) 
            && (Dew[ELEMENT(i+1,j-1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i+1,j-1)].GetWaterCourseValue()==noData)
            && (Dew[ELEMENT(i+1,j-1)].GetFlowLandValue()==128)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i+1,j-1)]);
          nUp++;
        }
        if ((j>0) 
            && (Dew[ELEMENT(i,j-1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i,j-1)].GetWaterCourseValue()==noData)
            && (Dew[ELEMENT(i,j-1)].GetFlowLandValue()==1)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i,j-1)]);
          nUp++;
        }
        if ((i>0 && j>0) 
            && (Dew[ELEMENT(i-1,j-1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i-1,j-1)].GetWaterCourseValue()==noData)
            && (Dew[ELEMENT(i-1,j-1)].GetFlowLandValue()==2)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i-1,j-1)]);
          nUp++;
        }
        if ((i>0)
            && (Dew[ELEMENT(i-1,j)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i-1,j)].GetWaterCourseValue()==noData)
            && (Dew[ELEMENT(i-1,j)].GetFlowLandValue()==4)) {
          Dew[ELEMENT(i,j)].SetUpLandFlow(nUp,&Dew[ELEMENT(i-1,j)]);
          nUp++;
        }
        if ((i>0 && j<nCols-1)
            && (Dew[ELEMENT(i-1,j+1)].GetLandIndex()!=noData)
            && (Dew[ELEMENT(i-1,j+1)].GetWaterCourseValue()==noData)
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


void WriteLandScapeAndWaterCourseFlow(DistributedElement * const Dew, int nRows, int nCols, int noData, ofstream &fout)
{
  int i,j;

  // Pointers between landscape elements  
  ofstream landUpFlow("dew_landupflow.txt");  // Open for writing
  for (i=0; i<nRows*nCols; i++) {
    if (Dew[i].GetNumUpLand()>0) {
      landUpFlow << Dew[i].GetLandIndex() << "  " << Dew[i].GetNumUpLand() << " :  ";
      for (j=0; j<Dew[i].GetNumUpLand(); j++) 
        landUpFlow << Dew[i].GetUpLandFlow(j)->GetLandIndex() << "  ";
      landUpFlow << endl;
    }
  }
  landUpFlow << endl;
  landUpFlow.close();
}


void DistanceSort(DistributedElement * const distElement, ParametersCommon * const ParCommonStore, 
                  MeteorologicalStations * const MetStations,
                  int * precipitationStation, int * temperatureStation, double * precipitationWeight, double * temperatureWeight)
{
  int i, j;
  int numP, minIndex;
  int temporaryStation;
  double temporaryDistance;
  double sumWeight;
  int * precipitationSort = new int [MetStations->GetNumPrecStations()];
  int * temperatureSort = new int [MetStations->GetNumTempStations()];
  double * distancePrec = new double [MetStations->GetNumPrecStations()];
  double * distanceTemp = new double [MetStations->GetNumTempStations()];

  for (i=0; i<MetStations->GetNumPrecStations(); i++) {
    precipitationSort[i] = i;
    distancePrec[i] = sqrt(pow((distElement->GetXCoord() - MetStations->GetStationCoordX(i)),2.0) 
                           + pow((distElement->GetYCoord() - MetStations->GetStationCoordY(i)),2.0));
    //    cout << "Prec " << precipitationSort[i] << "  " << distancePrec[i] << endl;
  }
  for (i=0; i<MetStations->GetNumPrecStations(); i++) {
    minIndex = i;
    for (j=i+1; j<MetStations->GetNumPrecStations(); j++) {
      if (distancePrec[j] < distancePrec[minIndex]) {
        temporaryDistance = distancePrec[minIndex];
        distancePrec[minIndex] = distancePrec[j];
        distancePrec[j] = temporaryDistance;
        temporaryStation = precipitationSort[minIndex];
        precipitationSort[minIndex] = precipitationSort[j];
        precipitationSort[j] = temporaryStation;
      }
    }
  }
  //  cout << "PREC  " << precipitationSort[0] << "  " << distancePrec[0] << endl;

  numP = MetStations->GetNumPrecStations();
  for (i=0; i<MetStations->GetNumTempStations(); i++) {
    temperatureSort[i] = i;
    distanceTemp[i] = sqrt(pow((distElement->GetXCoord() - MetStations->GetStationCoordX(numP+i)),2.0) 
                    + pow((distElement->GetYCoord() - MetStations->GetStationCoordY(numP+i)),2.0));
    //    cout << "Temp " << temperatureSort[i] << "  " << distanceTemp[i] << endl;
  }
  for (i=0; i<MetStations->GetNumTempStations(); i++) {
    minIndex = i;
    for (j=i+1; j<MetStations->GetNumTempStations(); j++) {
      if (distanceTemp[j] < distanceTemp[minIndex]) {
        temporaryDistance = distanceTemp[minIndex];
        distanceTemp[minIndex] = distanceTemp[j];
        distanceTemp[j] = temporaryDistance;
        temporaryStation = temperatureSort[minIndex];
        temperatureSort[minIndex] = temperatureSort[j];
        temperatureSort[j] = temporaryStation;
      }
    }
  }
  //  cout << "TEMP  "  << temperatureSort[0] << "  " << distanceTemp[0] << endl;

  sumWeight = 0.0;
  for (i=0; i<ParCommonStore->GetNUM_PREC_SERIES(); i++) {
    precipitationStation[i] = precipitationSort[i];
    precipitationWeight[i] = 1.0/distancePrec[i];
    sumWeight = sumWeight + precipitationWeight[i];
  }
  if (sumWeight != 1.0) {
    for (i=0; i<ParCommonStore->GetNUM_PREC_SERIES(); i++) {
      precipitationWeight[i] = precipitationWeight[i]/sumWeight;
    }
  }

  sumWeight = 0.0;
  for (i=0; i<ParCommonStore->GetNUM_TEMP_SERIES(); i++) {
    temperatureStation[i] = temperatureSort[i];
    temperatureWeight[i] = 1.0/distanceTemp[i];
    sumWeight = sumWeight + temperatureWeight[i];
  }
  if (sumWeight != 1.0) {
    for (i=0; i<ParCommonStore->GetNUM_TEMP_SERIES(); i++) {
      temperatureWeight[i] = temperatureWeight[i]/sumWeight;
    }
  }

  //  cout << distElement->GetLandIndex() << "  " << distElement->GetGeoIndex() << "  "
  //       << distElement->GetXCoord() << "  " << distElement->GetYCoord() << endl;

  delete [] precipitationSort;
  delete [] temperatureSort;
  delete [] distancePrec;
  delete [] distanceTemp;

}
