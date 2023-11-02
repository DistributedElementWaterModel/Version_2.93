/****************************************************************************************
 *                                                                                      *
 *  Program for calculating water balance for distributed hydrological response units,  *
 *  and traversing and accumulating runoff from a hierarchy of sub-catchments,          *
 *  river networks, lakes and landscape elements.                                       *
 *                                                                                      *
 *  Author: Stein Beldring, Norway                                                      *
 *                                                                                      *
 ****************************************************************************************/

#include "date_time.h"
#include "classDew.h"

void WaterBalanceTimeSeries(DistributedElement * const Dew, ParametersCommon * const ParCommonStore,
                            MeteorologicalStations * const MetStations, DateTimeInformation * const DateTimeStore,
                            InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore, 
                            int numLand, int timeStep, int initialTimeSteps, bool * inputDataFound, bool * firstTotal);
void WaterBalanceGrid(DistributedElement * const Dew,  ParametersCommon * ParCommonStore, 
                      InputElement * InputElementStore, DateTimeInformation * const DateTimeStore, 
                      int numLand, int timeStep, int initialTimeSteps, int nRows, int nCols, DateTime datetime, char * metPath,
                      unsigned short int * precip10, unsigned short int * temp10K, bool flowHierarchy, 
                      bool * inputDataFound, bool * firstTotal);
void SnowGlacierIceReDistribution(WaterCourse ** const Outlet, DistributedElement * const Dew, ParametersCommon * ParCommonStore, int initialTimeSteps, int numberTimeSteps,
                                  int numLand, int numWatcOut, int timeStep, DateTime datetime, 
                                  int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, bool modelCalibration, ofstream &fout);
void SetDistributedElementInput(DistributedElement * const thisElement, ParametersCommon * const ParCommonStore,
                            MeteorologicalStations * const MetStations, 
                            InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore, 
                            int timeStep, bool * inputDataFound);
void ReadWaterCourseIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, 
                               int numWatc, bool flowHierarchy, ifstream &fileControl, ofstream &fout);
void ReadLandscapeHierarchy(DistributedElement * const Dew, ifstream &fileControl, ofstream &fout);
void WriteWaterCourseIdentifier(WaterCourse * const WaterElement, int numWatc, ofstream &fout);
void TraverseCorrectionWaterCourse(WaterCourse * const thisWaterCourse, int numberCorrectionCatchments,
                                   int * correctionCatchments, double * correctionPrecipitation, 
                                   double * correctionTemperature);
void TraverseCorrectionLandScape(DistributedElement * const thisElement, double precCorr, double tempCorr);
// ** Function to be performed in case: no input to landscape element from upstream elements
/*void TraverseWaterCourse(WaterCourse * const thisWaterCourse, int timeStep, ofstream &fout);
  void TraverseLandScape(DistributedElement * const thisElement, int timeStep, ofstream &fout);*/
// ** Function to be performed in case: input to landscape element from upstream elements
void TraverseWaterCourse(WaterCourse * const thisWaterCourse, ParametersCommon * const ParCommonStore, 
  MeteorologicalStations * const MetStations, InputTimeSeries * InputTimeSeriesStore,
  InputElement * InputElementStore, DateTimeInformation * const DateTimeStore, int timeStep, int initialTimeSteps, char inputFormat,
  bool flowHierarchy, bool forceDirect, bool * firstTotal, ofstream &fout);
void TraverseLandScape(DistributedElement * const thisElement, ParametersCommon * const ParCommonStore, 
  MeteorologicalStations * const MetStations, InputTimeSeries * InputTimeSeriesStore,
  InputElement * InputElementStore, DateTimeInformation * const DateTimeStore, int timeStep, int initialTimeSteps, char inputFormat,
  bool * waterCourseInputFound, bool flowHierarchy, bool forceDirect, bool * firstTotal, ofstream &fout);
// ** End algorithms to be performed in different cases of input to landscape element
void TraverseMissingDataWaterCourse(WaterCourse * const thisWaterCourse, int timeStep, ofstream &fout);
void SetMissingDataWaterCourse(WaterCourse * const thisWaterCourse, int timeStep, ofstream &fout);
void SetMissingDataLandScape(DistributedElement * const thisElement, int timeStep, ofstream &fout);
void WriteWaterCourseDischarge(WaterCourse * WaterElement, int numWatc, DateTimeInformation * const DateTimeStore,
                                int secondsPerTimeStep, bool modelCalibration, ofstream &fout);
void WriteWaterCourseWaterBalance(WaterCourse * WaterElement, int numWatc, DateTimeInformation * const DateTimeStore, 
                                   int secondsPerTimeStep);
void WriteDistributedElementTimeSeries(DistributedElement * const Dew, int numLand, DateTimeInformation * const DateTimeStore,
                                       int secondsPerTimeStep);
void WriteOutletDischarge(WaterCourse * const thisWaterCourse, DateTime startSimulationTime, 
                          DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                          int secondsPerTimeStep, bool modelCalibration, ofstream &fout);
void WriteOutletWaterBalance(WaterCourse * const thisWaterCourse, DateTime startSimulationTime, 
                             DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                             int secondsPerTimeStep, bool flowHierarchy, bool forceDirect);
void WriteTotalReservoirStorage(TotalReservoirStorage * const TotalReservoirStore, DateTimeInformation * const DateTimeStore,
                                bool modelCalibration, int secondsPerTimeStep);
void ObjectiveCriteria(int numberTimeSteps, double *obs, double *sim, 
                       double *ns, double *rmse, double *bias, double *pears);
void WriteReducedBinaryGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,int nCols, 
                    int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void WriteBinaryGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,int nCols, 
                    int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void WriteAsciiGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,int nCols, 
                    int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout);
void SetCommonParameters(ParametersCommon * const ParCommonStore, ifstream &fileControl, ofstream &fout);
void SetLandSurfaceParameters(ParametersLandSurface * const ParLandSurfaceStore, ifstream &fileControl, ofstream &fout);
void SetSnowDistribution(ParametersLandSurface * const thisParLandSurface, double cvSnow);
void SetSubSurfaceHbvParameters(ParametersSubSurfaceHbv * const ParSubSurfaceHbvStore, ifstream &fileControl, ofstream &fout);
void SetKiWaParameters(ParametersKiWa * const ParKiWaStore, ifstream &fileControl, ofstream &fout);

double power(double, double);
int leapYear(int year);
int dayNumber(int year, int month, int day);
void dayNo2Date(int dayNo, int year, int * month, int * day);
double potentialEvap (double, double);
double HBVTranspSoilEvap(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel);
double KiWaTranspSoilEvap(double teta, double temp, double epotPar, double eactPar, double tSat0, double wiltPoint);

int main(int argc, char *argv[])
{
  std::cout << "\n\n Distributed Element Water Model \n\n";

  bool inputDataFound=false;
  bool modelCalibration=false;
  bool flowHierarchy=false;
  bool forceDirect=false;
  bool firstTotal=false;
  bool areaFractionError=false;
  bool lakeAreaFound=false;
  bool glacierAreaFound=false;
  bool landAreaFound=false;
  char fileName[100];
  char fileNameInput[100];
  char fileNameWaterCourse[100];
  char buffer[256];
  char modelRun = 'R';
  char hierarchy = 'H';
  char inputFormat = 'F';
  char ch;
  char *metPath;
  int i, j, k;
  int numLand;
  int geoIndex;
  int landIndex;
  int nRows, nCols, noData;
  int modStruct;
  int landSurf;
  int soil;
  int elementFlowDir;
  int timeStep;
  int numWatc, numWatcUp, numWatcOut;
  int waterCourseId;
  int numberCorrectionCatchments;
  int seriesNumber;
  //  int time, startModelTime, startSimulationTime, endSimulationTime;
  int initialTimeSteps;
  int numberTimeSteps;
  int lakeNumber;
  double seriesWeight;
  double sumWeight;
  double sumArea;
  double precStationsWeightedElevation, tempStationsWeightedElevation;
  double correction;
  double xllCorner, yllCorner, cellSize;
  double elementArea, elementElevation, elementSlopeLength, elementSlopeAngle, elementAspect, lakePercent, glacierPercent;
  double areaFraction[maximumNumberLandClasses];
  DateTime startModelTime;
  DateTime startSimulationTime;
  DateTime endSimulationTime;
  DateTime datetime;
  Lake *ptrLake;
  Glacier *ptrGlacier;
  HbvAquifer *ptrHbvAquifer;
  HbvAquifer *lastHbvAquifer;
  KwaAquifer *ptrKwaAquifer;
  KwaAquifer *lastKwaAquifer;
  LakeWaterBalance *ptrLakeWaterBalance;
  Vegetation *ptrVegetation;
  Snow *ptrSnow;
  GlacierSurface * ptrSurface;
  //  GlacierIce *ptrIce;
  HBV *ptrHbv;
  KinematicWave *ptrKwa;
  MODEL_STRUCTURE modelStructure;
  LANDSURFACE landSurfType[numberLandSurfaceClasses];
  SOIL soilType[numberLandSurfaceClasses];

  DateTime nowDate;
  nowDate.now();
  DateTime finalDate(finalYear,finalMonth,finalDay,finalHour,finalMinute,0);
  //  cout<<" "<<finalDate.getYear()<<" "<<finalDate.getMonth()<<" "<<finalDate.getDay()<<" "<<finalDate.getHour()<<" "<<finalDate.getMinute()<<endl;
/*  if (nowDate > finalDate) {
    strcpy(buffer,"del ");
    strcat(buffer,argv[0]);
    strcat(buffer,"*\n");
    strcat(buffer,"del dew.bat\n");
    ofstream fexit("abc_987.bat");
    fexit << "@echo off\n" << buffer << endl;
    fexit.close();
    system("copy abc_987.bat dew.bat /Y >NUL");
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
 
  /*  while (modelRun != 'S' && modelRun != 'C' && modelRun != 's' && modelRun != 'c') {
      cout << " Type of model run, simulation(S) or calibration(C): ";
      cin >> modelRun;
      cout << endl;
      }*/   
  fileControl.ignore(100,':');
  fileControl >> modelRun;
  fileControl.ignore(256,'\n');
  if (modelRun != 'S' && modelRun != 'C' && modelRun != 's' && modelRun != 'c') {
    cout << "\n Type of model run, simulation(S) or calibration(C) \n\n";
    exit (1);
  }
  if (modelRun == 'C' || modelRun == 'c') modelCalibration = true;

  /*  while (hierarchy != 'N' && hierarchy != 'C' && hierarchy != 'n' && hierarchy != 'c') {
      cout << " Landscape elements hierarchy, flow direction network(N) or nested catchments(C): ";
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
  if (hierarchy == 'N' || hierarchy == 'n') flowHierarchy = true;

  /*  while (inputFormat != 'G' && inputFormat != 'T' && inputFormat != 'g' && inputFormat != 't') {
      cout << " Input data format, grid files(G) or time series file(T): ";
      cin >> inputFormat;
      cout << endl;
      }*/   
  fileControl.ignore(100,':');
  fileControl >> inputFormat;
  fileControl.ignore(256,'\n');
  if (inputFormat != 'G' && inputFormat != 'T' && inputFormat != 'g' && inputFormat != 't') {
    cout << "\n Input data format, grid files(G) or time series file(T) \n\n";
    exit (1);
  }

  /*  cout << " Output file: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ofstream fout(fileName);  // Open for writing

  // Read time steps for start model, start simulation, end simulation
  /*  cout << " Start model time : ";
      cin >> startModelTime;
      cout << endl;
      cout << " Start simulation time: ";
      cin >> startSimulationTime;
      cout << endl;
      cout << " End simulation time: ";
      cin >> endSimulationTime;
      cout << endl;
      if (startModelTime > startSimulationTime || startSimulationTime > endSimulationTime) {
      cout << endl << " Time step error " << endl << endl;
      exit(1);
      }
      initialTimeSteps = (int)(startSimulationTime-startModelTime);
      numberTimeSteps = (int)(endSimulationTime-startSimulationTime+1);
      cout << " initialTimeSteps = " << initialTimeSteps << " numberTimeSteps = "
      << numberTimeSteps << endl;*/

  // Object for storing date and time information
  DateTimeInformation * DateTimeStore = new DateTimeInformation;
  DateTimeStore->ReadDateTimeInformation(fileControl);

  // Object for storing meteorological station information
  MeteorologicalStations * MetStations = new MeteorologicalStations;
  MetStations->SetMeteorologicalStations(fileControl, fout);

  // Read common parameters file and set parameter values
  ParametersCommon * ParCommonStore = new ParametersCommon;
  SetCommonParameters(ParCommonStore, fileControl, fout);
  if ((ParCommonStore->GetSECONDS_TIMESTEP() % minimumTimeStep) != 0) {
    cout << endl << " Temporal resolution must be a multiple of " << minimumTimeStep << " seconds : " << ParCommonStore->GetSECONDS_TIMESTEP() << endl << endl;
    exit(1);
  }
  // End read common parameters

  // Store data in object for date and time information: start model, start simulation, end simulation,
  // no. of initial time steps and no. simulation time steps
  DateTimeStore->SetCommonPar(ParCommonStore);
  DateTimeStore->CalculateTimeStepInformation();
  DateTimeStore->WriteDateTimeInformation(fout);
  initialTimeSteps = DateTimeStore->GetInitialTimeSteps();
  numberTimeSteps = DateTimeStore->GetNumberTimeSteps();
  startModelTime = DateTimeStore->GetStartModelTime();
  startSimulationTime = DateTimeStore->GetStartSimulationTime();
  endSimulationTime = DateTimeStore->GetEndSimulationTime();

  // Read landsurface parameters file and set parameter values
  ParametersLandSurface * ParLandSurfaceStore = new ParametersLandSurface [numberLandSurfaceClasses];
  SetLandSurfaceParameters(ParLandSurfaceStore, fileControl, fout);
  // End read landsurface parameters

  // Read subsurface parameters file and set parameter values
  ParametersSubSurfaceHbv * ParSubSurfaceHbvStore = new ParametersSubSurfaceHbv [numberSoilClasses];
  SetSubSurfaceHbvParameters(ParSubSurfaceHbvStore, fileControl, fout);
  // End read subsurface parameters

  // Read kinematic wave parameters file and set parameter values
  ParametersKiWa * ParKiWaStore = new ParametersKiWa [numberSoilClasses];
  SetKiWaParameters(ParKiWaStore, fileControl, fout);
  // End read kinematic wave parameters

  // Object for storing landscape elements selected for HBV time series output
  SelectedHbvTimeSeriesElements * SelectedHbvTimeSeriesElementsStore = new SelectedHbvTimeSeriesElements;
  SelectedHbvTimeSeriesElementsStore->SelectedHbvTimeSeriesElementsInput(fileControl, fout);

  // Object for storing landscape elements selected for kinematic wave time series output
  SelectedKiWaTimeSeriesElements * SelectedKiWaTimeSeriesElementsStore = new SelectedKiWaTimeSeriesElements;
  SelectedKiWaTimeSeriesElementsStore->SelectedKiWaTimeSeriesElementsInput(fileControl, fout);

  // Object for storing landscape elements selected for kinematic wave hillslope output
  SelectedKiWaHillslopeElements * SelectedKiWaHillslopeElementsStore = new SelectedKiWaHillslopeElements;
  SelectedKiWaHillslopeElementsStore->SelectedKiWaHillslopeElementsInput(fileControl, fout);

  // Object for storing watercourse/catchment elements selected for time series output
  SelectedWaterCourseTimeSeriesElements * SelectedWaterCourseTimeSeriesElementsStore = new SelectedWaterCourseTimeSeriesElements;
  SelectedWaterCourseTimeSeriesElementsStore->SelectedWaterCourseTimeSeriesElementsInput(fileControl, fout);

  // Object for storing input data for each landscape element
  InputElement * InputElementStore = new InputElement(numberInputSeries);

  // Define total reservoir storage object
  TotalReservoirStorage * TotalReservoirStore = new TotalReservoirStorage;
  TotalReservoirStore->AllocateTotalReservoirStorage(initialTimeSteps+numberTimeSteps);
  TotalReservoirStore->SetCommonPar(ParCommonStore);
  TotalReservoirStore->SetInitialTotalReservoirStorage();
  // End define total reservoir storage object

  // Read landscape element file and generate landscape element objects  
  /*  cout << " File with landscape element information: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finLandScape(fileName);  // Open for reading
  if (finLandScape == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  } 
  finLandScape >> buffer >> nCols;
  finLandScape >> buffer >> nRows;
  finLandScape >> buffer >> xllCorner;
  finLandScape >> buffer >> yllCorner;
  finLandScape >> buffer >> cellSize;
  finLandScape >> buffer >> noData;
  finLandScape.ignore(100,':');
  finLandScape >> numLand;
  DistributedElement * Dew = new DistributedElement [numLand];
  if (!Dew) {
    cout << endl << " Error during dynamic memory allocating for Dew[" << numLand << "]" << endl << endl;
    exit(1);
  } 

  for (i=0; i<numLand; i++) {
    finLandScape >> landIndex >> geoIndex >> modStruct >> elementArea >> elementElevation;
    finLandScape >> elementSlopeLength >> elementSlopeAngle >> elementAspect >> elementFlowDir >> lakePercent >> glacierPercent;
    modelStructure=MODEL_STRUCTURE(modStruct); 
    //    cout << landIndex << "  " << geoIndex << "  " << elementArea << "  " << elementElevation << "  " << elementSlopeLength << "  " << elementSlopeAngle << "  " << elementAspect << "  " << elementFlowDir << "  " << lakePercent << "  " << glacierPercent << "  " << endl;
    sumArea = lakePercent + glacierPercent;
    for (j=0; j<maximumNumberLandClasses; j++) {
      finLandScape >> landSurf >> soil >> areaFraction[j]; 
      landSurfType[j]=LANDSURFACE(landSurf);
      soilType[j]=SOIL(soil);
      sumArea = sumArea + areaFraction[j];
    }
    //for (j=0; j<maximumNumberLandClasses; j++) {
    //  cout << landSurfType[j] << "  " << soilType[j] << "  " << areaFraction[j] << "  "; 
    //}
    //cout << endl;
    // Explicit flow hierarchy: control number of land cover types in each landscape element
    if (flowHierarchy) {
      areaFractionError=false;
      if (lakePercent > 0.0) lakeAreaFound = true;
      else lakeAreaFound = false;    
      if (glacierPercent > 0.0) glacierAreaFound = true;
      else glacierAreaFound = false;    
      landAreaFound = false;
      for (j=0; j<maximumNumberLandClasses; j++) {
        if (areaFraction[j] > 0.0) landAreaFound = true;
      }
      if (lakePercent == 100.0) {
//      cout << "lakePercent\n";
        if (glacierPercent != 0.0) areaFractionError = true;
        for (j=0; j<maximumNumberLandClasses; j++) {
          if (areaFraction[j] != 0.0) areaFractionError = true; 
        }
      }
      if (glacierPercent == 100.0) {
//      cout << "glacierPercent\n";
        if (lakePercent != 0.0) areaFractionError = true;
        for (j=0; j<maximumNumberLandClasses; j++) {
          if (areaFraction[j] != 0.0) areaFractionError = true; 
        }
      }
      for (j=0; j<maximumNumberLandClasses; j++) {
        if (areaFraction[j] == 100.0) {
//        cout << "areaFraction[" << j << "]" << endl;
          if (lakePercent != 0.0) areaFractionError = true;
          if (glacierPercent != 0.0) areaFractionError = true;
          for (k=0; k<maximumNumberLandClasses; k++) {
            if (k != j && areaFraction[k] != 0.0) areaFractionError = true; 
          }
        }
      }
//      cout << areaFractionError << "  " << lakeAreaFound << "  " << glacierAreaFound << "  " << landAreaFound << endl;
      if (lakeAreaFound && (glacierAreaFound || landAreaFound)) areaFractionError = true;
      if (glacierAreaFound && (lakeAreaFound || landAreaFound)) areaFractionError = true;
      if (landAreaFound && (lakeAreaFound || glacierAreaFound)) areaFractionError = true;
      if (areaFractionError) {
        cout << endl << " One land cover type in each landscape element for explicit flow hierarchy: land index ";
        cout << landIndex << endl << endl;
        exit (1);
      } 
    }
    // End control number of land cover types in each landscape element
    if (sumArea != 100.0) {
      lakePercent = lakePercent*100.0/sumArea;
      glacierPercent = glacierPercent*100.0/sumArea;
      for (j=0; j<maximumNumberLandClasses; j++) 
        areaFraction[j] = areaFraction[j]*100.0/sumArea;
    }
    //cout << lakePercent << "  " << glacierPercent << "  ";
    //for (j=0; j<maximumNumberLandClasses; j++) {
    //    cout << areaFraction[j] << "  "; 
    //}
    //cout << endl;
    Dew[landIndex].SetGeoIndex(geoIndex);
    Dew[landIndex].SetLandIndex(landIndex);
    Dew[landIndex].SetArea(elementArea);
    Dew[landIndex].SetElevation(elementElevation);
    Dew[landIndex].SetSlopeLength(elementSlopeLength);
    Dew[landIndex].SetSlopeAngle(elementSlopeAngle);
    Dew[landIndex].SetAspect(elementAspect);
    Dew[landIndex].SetFlowDirection(elementFlowDir);
    Dew[landIndex].SetSelectedKiWaTimeSeriesElements(SelectedKiWaTimeSeriesElementsStore);
    Dew[landIndex].SetSelectedHbvTimeSeriesElements(SelectedHbvTimeSeriesElementsStore);
    Dew[landIndex].SetCommonPar(ParCommonStore);
    Dew[landIndex].SetTotalReservoir(TotalReservoirStore);
    Dew[landIndex].AllocateInputValues(InputElementStore->GetNumberValues());
    // Allocate space for water balance time series for landscape elements 
    for (j=0; j<Dew[landIndex].GetSelectedKiWaTimeSeriesElements()->GetNumberElements(); j++) {
      if (Dew[landIndex].GetLandIndex() == Dew[landIndex].GetSelectedKiWaTimeSeriesElements()->GetKiWaTimeSeriesElement(j)) {
        Dew[landIndex].AllocateKiWaWaterBalance(initialTimeSteps+numberTimeSteps);
      }
    }
    for (j=0; j<Dew[landIndex].GetSelectedHbvTimeSeriesElements()->GetNumberElements(); j++) {
      if (Dew[landIndex].GetLandIndex() == Dew[landIndex].GetSelectedHbvTimeSeriesElements()->GetHbvTimeSeriesElement(j)) {
        Dew[landIndex].AllocateHbvWaterBalance(initialTimeSteps+numberTimeSteps);
      }
    }
    // Time series input data format - read meteorological station information
    if (inputFormat == 'T' || inputFormat == 't') {
      sumWeight=0.0;
      precStationsWeightedElevation=0.0;
      Dew[landIndex].AllocateMetSeries(ParCommonStore->GetNUM_PREC_SERIES(),ParCommonStore->GetNUM_TEMP_SERIES());
      for (j=0; j<ParCommonStore->GetNUM_PREC_SERIES(); j++) {
        finLandScape >> seriesNumber >> seriesWeight; 
        Dew[landIndex].SetMetSeriesNumber(j, seriesNumber);
        Dew[landIndex].SetMetSeriesWeight(j, seriesWeight);
        precStationsWeightedElevation = precStationsWeightedElevation+MetStations->GetStationAltitude(seriesNumber)*seriesWeight;
        sumWeight = sumWeight+seriesWeight;
      }
      if (sumWeight != 1.0) {
        for (j=0; j<ParCommonStore->GetNUM_PREC_SERIES(); j++) {
          seriesWeight = Dew[landIndex].GetMetSeriesWeight(j)/sumWeight;
          Dew[landIndex].SetMetSeriesWeight(j, seriesWeight);
          precStationsWeightedElevation = precStationsWeightedElevation/sumWeight;
        }
      }
      //      cout << " precStationsWeightedElevation = " << precStationsWeightedElevation << endl;
      sumWeight=0.0;
      tempStationsWeightedElevation=0.0;
      for (j=0; j<ParCommonStore->GetNUM_TEMP_SERIES(); j++) {
        finLandScape >> seriesNumber >> seriesWeight; 
        Dew[landIndex].SetMetSeriesNumber(ParCommonStore->GetNUM_PREC_SERIES()+j,seriesNumber);
        Dew[landIndex].SetMetSeriesWeight(ParCommonStore->GetNUM_PREC_SERIES()+j, seriesWeight);
        tempStationsWeightedElevation = tempStationsWeightedElevation+MetStations->GetStationAltitude(seriesNumber)*seriesWeight;
        sumWeight = sumWeight+seriesWeight;
      }
      if (sumWeight != 1.0) {
        for (j=0; j<ParCommonStore->GetNUM_TEMP_SERIES(); j++) {
          seriesWeight = Dew[landIndex].GetMetSeriesWeight(ParCommonStore->GetNUM_PREC_SERIES()+j)/sumWeight;
          Dew[landIndex].SetMetSeriesWeight(ParCommonStore->GetNUM_PREC_SERIES()+j, seriesWeight);
          tempStationsWeightedElevation = tempStationsWeightedElevation/sumWeight;
        }
      }
      //      cout << " tempStationsWeightedElevation = " << tempStationsWeightedElevation << endl;
        
      /*      for (j=0; j<ParCommonStore->GetNUM_PREC_SERIES(); j++) {
              cout << "P " << Dew[landIndex].GetMetSeriesNumber(j) << "  " << Dew[landIndex].GetMetSeriesWeight(j) << endl;
              }
              for (j=0; j<ParCommonStore->GetNUM_TEMP_SERIES(); j++) {
              cout << "T " << Dew[landIndex].GetMetSeriesNumber(ParCommonStore->GetNUM_PREC_SERIES()+j) << 
              "  " << Dew[landIndex].GetMetSeriesWeight(ParCommonStore->GetNUM_PREC_SERIES()+j) << endl;
              }*/
    }
    else {    
      finLandScape.ignore(256,'\n');
    }
    if (lakePercent>0.0) {
      ptrLake = new Lake;
      ptrLake->SetAreaFraction(lakePercent);
      ptrLakeWaterBalance = new LakeWaterBalance;
      ptrLakeWaterBalance->SetCommonPar(ParCommonStore);
      //      ptrLakeWaterBalance->SetInputTimeSeries(InputTimeSeriesStore);
      ptrLakeWaterBalance->SetInputElement(InputElementStore);
      ptrLakeWaterBalance->SetLandScapeElement(&Dew[landIndex]);
      ptrLakeWaterBalance->SetLakeValues(ParCommonStore->GetINITIAL_LAKE_TEMP(), ParCommonStore->GetINITIAL_LAKE_LEVEL());
      ptrLake->SetLakeWaterBalance(ptrLakeWaterBalance);
      Dew[landIndex].SetLake(ptrLake);
    }
    if (glacierPercent>0.0) {
      ptrGlacier = new Glacier;
      ptrGlacier->SetAreaFraction(glacierPercent);
      ptrSurface = new GlacierSurface;
      ptrSurface->SetCommonPar(ParCommonStore);
      ptrSurface->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]); 
      //      ptrSurface->SetInputTimeSeries(InputTimeSeriesStore);
      ptrSurface->SetInputElement(InputElementStore);
      ptrSurface->SetLandScapeElement(&Dew[landIndex]);
      ptrGlacier->SetGlacierSurface(ptrSurface);
      ptrSnow = new Snow;
      ptrSnow->SetCommonPar(ParCommonStore);
      ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]); 
      //      ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
      ptrSnow->SetInputElement(InputElementStore);
      ptrSnow->SetLandScapeElement(&Dew[landIndex]);
      ptrSnow->SetSnowStore(ParCommonStore->GetINITIAL_SNOW());
      ptrGlacier->SetSnow(ptrSnow);
      //      ptrIce = new GlacierIce;
      //      ptrIce->SetCommonPar(ParCommonStore);
      //      ptrIce->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]); 
      //      ptrIce->SetInputTimeSeries(InputTimeSeriesStore);
      //      ptrIce->SetInputElement(InputElementStore);
      //      ptrIce->SetLandScapeElement(&Dew[landIndex]);
      //      ptrGlacier->SetGlacierIce(ptrIce);
      if (modelStructure==HBV_MODEL) {
        ptrHbv = new HBV;
        ptrHbv->SetCommonPar(ParCommonStore);
        ptrHbv->SetSubSurfaceHbvPar(&ParSubSurfaceHbvStore[GLACIER_BED]);
        ptrHbv->SetLandSurfacePar(&ParLandSurfaceStore[GLACIER]); 
        //      ptrHbv->SetInputTimeSeries(InputTimeSeriesStore);
        ptrHbv->SetInputElement(InputElementStore);
        ptrHbv->SetLandScapeElement(&Dew[landIndex]);
        ptrHbv->SetInitialHbvValues();
        ptrHbv->SetSubSurfaceHbvStore(ParCommonStore->GetINITIAL_SOIL_MOISTURE(), 
                                      ParCommonStore->GetINITIAL_UPPER_ZONE(), ParCommonStore->GetINITIAL_LOWER_ZONE());
        ptrGlacier->SetHBV(ptrHbv);
      }
      else if (modelStructure==KWA_MODEL) { 
        ptrKwa = new KinematicWave;
        ptrKwa->SetCommonPar(ParCommonStore);
        ptrKwa->SetKiWaPar(&ParKiWaStore[soilType[GLACIER_BED]]);
        ptrKwa->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[GLACIER]]); 
        //      ptrKwa->SetInputTimeSeries(InputTimeSeriesStore);
        ptrKwa->SetInputElement(InputElementStore);
        ptrKwa->SetLandScapeElement(&Dew[landIndex]);
        ptrKwa->SetInitialKwaValues();
        ptrGlacier->SetKinematicWave(ptrKwa);
      }       
      Dew[landIndex].SetGlacier(ptrGlacier);
    }
    if (modelStructure==HBV_MODEL) {
      if (areaFraction[0]>0.0) {
        ptrHbvAquifer = new HbvAquifer;
        ptrHbvAquifer->SetAreaFraction(areaFraction[0]);
        ptrVegetation = new Vegetation;
        ptrVegetation->SetCommonPar(ParCommonStore);
        ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]); 
        //      ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
        ptrVegetation->SetInputElement(InputElementStore);
        ptrVegetation->SetLandScapeElement(&Dew[landIndex]);
        ptrHbvAquifer->SetVegetation(ptrVegetation);
        ptrSnow = new Snow;
        ptrSnow->SetCommonPar(ParCommonStore);
        ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]); 
        //      ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
        ptrSnow->SetInputElement(InputElementStore);
        ptrSnow->SetLandScapeElement(&Dew[landIndex]);
        ptrSnow->SetSnowStore(ParCommonStore->GetINITIAL_SNOW());
        ptrHbvAquifer->SetSnow(ptrSnow);
        ptrHbv = new HBV;
        ptrHbv->SetCommonPar(ParCommonStore);
        ptrHbv->SetSubSurfaceHbvPar(&ParSubSurfaceHbvStore[soilType[0]]);
        ptrHbv->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]); 
        //      ptrHbv->SetInputTimeSeries(InputTimeSeriesStore);
        ptrHbv->SetInputElement(InputElementStore);
        ptrHbv->SetLandScapeElement(&Dew[landIndex]);
        ptrHbv->SetInitialHbvValues();
        ptrHbv->SetSubSurfaceHbvStore(ParCommonStore->GetINITIAL_SOIL_MOISTURE(), 
                                      ParCommonStore->GetINITIAL_UPPER_ZONE(), ParCommonStore->GetINITIAL_LOWER_ZONE());
        ptrHbvAquifer->SetHBV(ptrHbv);
        Dew[landIndex].SetHbvAquifer(ptrHbvAquifer);
        j=1;
        while (j<maximumNumberLandClasses && areaFraction[j]>0.0) {
          lastHbvAquifer=ptrHbvAquifer;
          ptrHbvAquifer = new HbvAquifer;
          ptrHbvAquifer->SetAreaFraction(areaFraction[j]);
          ptrVegetation = new Vegetation;
          ptrVegetation->SetCommonPar(ParCommonStore);
          ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]); 
          //        ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
          ptrVegetation->SetInputElement(InputElementStore);
          ptrVegetation->SetLandScapeElement(&Dew[landIndex]);
          ptrHbvAquifer->SetVegetation(ptrVegetation);
          ptrSnow = new Snow;
          ptrSnow->SetCommonPar(ParCommonStore);
          ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]); 
          //        ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
          ptrSnow->SetInputElement(InputElementStore);
          ptrSnow->SetLandScapeElement(&Dew[landIndex]);
          ptrSnow->SetSnowStore(ParCommonStore->GetINITIAL_SNOW());
          ptrHbvAquifer->SetSnow(ptrSnow);
          ptrHbv = new HBV;
          ptrHbv->SetCommonPar(ParCommonStore);
          ptrHbv->SetSubSurfaceHbvPar(&ParSubSurfaceHbvStore[soilType[j]]); 
          ptrHbv->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]); 
          //        ptrHbv->SetInputTimeSeries(InputTimeSeriesStore);
          ptrHbv->SetInputElement(InputElementStore);
          ptrHbv->SetLandScapeElement(&Dew[landIndex]);
          ptrHbv->SetInitialHbvValues();
          ptrHbv->SetSubSurfaceHbvStore(ParCommonStore->GetINITIAL_SOIL_MOISTURE(), 
                                        ParCommonStore->GetINITIAL_UPPER_ZONE(), ParCommonStore->GetINITIAL_LOWER_ZONE());
          ptrHbvAquifer->SetHBV(ptrHbv);
          lastHbvAquifer->SetNextHbvAquifer(ptrHbvAquifer);
          j++;
        }
      }
    }
    else if (modelStructure==KWA_MODEL) { 
      if (areaFraction[0]>0.0) {
        ptrKwaAquifer = new KwaAquifer;
        ptrKwaAquifer->SetAreaFraction(areaFraction[0]);
        ptrVegetation = new Vegetation;
        ptrVegetation->SetCommonPar(ParCommonStore);
        ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]); 
        //      ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
        ptrVegetation->SetInputElement(InputElementStore);
        ptrVegetation->SetLandScapeElement(&Dew[landIndex]);
        ptrKwaAquifer->SetVegetation(ptrVegetation);
        ptrSnow = new Snow;
        ptrSnow->SetCommonPar(ParCommonStore);
        ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]); 
        //      ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
        ptrSnow->SetInputElement(InputElementStore);
        ptrSnow->SetLandScapeElement(&Dew[landIndex]);
        ptrSnow->SetSnowStore(ParCommonStore->GetINITIAL_SNOW());
        ptrKwaAquifer->SetSnow(ptrSnow);
        ptrKwa = new KinematicWave;
        ptrKwa->SetCommonPar(ParCommonStore);
        ptrKwa->SetKiWaPar(&ParKiWaStore[soilType[0]]);
        ptrKwa->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[0]]); 
        //      ptrKwa->SetInputTimeSeries(InputTimeSeriesStore);
        ptrKwa->SetInputElement(InputElementStore);
        ptrKwa->SetDateTimeInfo(DateTimeStore);
        ptrKwa->SetSelectedKiWaHillslopeElements(SelectedKiWaHillslopeElementsStore);
        ptrKwa->SetLandScapeElement(&Dew[landIndex]);
        ptrKwa->SetInitialKwaValues();
        ptrKwaAquifer->SetKinematicWave(ptrKwa);
        Dew[landIndex].SetKwaAquifer(ptrKwaAquifer);
        j=1;
        while (j<maximumNumberLandClasses && areaFraction[j]>0.0) {
          lastKwaAquifer=ptrKwaAquifer;
          ptrKwaAquifer = new KwaAquifer;
          ptrKwaAquifer->SetAreaFraction(areaFraction[j]);
          ptrVegetation = new Vegetation;
          ptrVegetation->SetCommonPar(ParCommonStore);
          ptrVegetation->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]); 
          //        ptrVegetation->SetInputTimeSeries(InputTimeSeriesStore);
          ptrVegetation->SetInputElement(InputElementStore);
          ptrVegetation->SetLandScapeElement(&Dew[landIndex]);
          ptrKwaAquifer->SetVegetation(ptrVegetation);
          ptrSnow = new Snow;
          ptrSnow->SetCommonPar(ParCommonStore);
          ptrSnow->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]); 
          //        ptrSnow->SetInputTimeSeries(InputTimeSeriesStore);
          ptrSnow->SetInputElement(InputElementStore);
          ptrKwa->SetDateTimeInfo(DateTimeStore);
          ptrKwa->SetSelectedKiWaHillslopeElements(SelectedKiWaHillslopeElementsStore);
          ptrSnow->SetLandScapeElement(&Dew[landIndex]);
          ptrSnow->SetSnowStore(ParCommonStore->GetINITIAL_SNOW());
          ptrKwaAquifer->SetSnow(ptrSnow);
          ptrKwa = new KinematicWave;
          ptrKwa->SetCommonPar(ParCommonStore);
          ptrKwa->SetKiWaPar(&ParKiWaStore[soilType[j]]); 
          ptrKwa->SetLandSurfacePar(&ParLandSurfaceStore[landSurfType[j]]); 
          //        ptrKwa->SetInputTimeSeries(InputTimeSeriesStore);
          ptrKwa->SetInputElement(InputElementStore);
          ptrKwa->SetLandScapeElement(&Dew[landIndex]);
          ptrKwa->SetInitialKwaValues();
          ptrKwaAquifer->SetKinematicWave(ptrKwa);
          lastKwaAquifer->SetNextKwaAquifer(ptrKwaAquifer);
          j++;
        }
      }
    }
  }

  ptrLake=0;
  ptrGlacier=0;
  ptrHbvAquifer=0;
  lastHbvAquifer=0;
  ptrKwaAquifer=0;
  lastKwaAquifer=0;
  ptrVegetation=0;
  ptrSnow=0;
  //  ptrIce=0;
  ptrHbv=0;
  ptrKwa=0;
  finLandScape.close();
  // End read file with landscape elements 

  // Read watercourse information file and generate watercourse objects  
  /*  cout << " File with watercourse/catchment hierarchy: ";
      cin >> fileNameWaterCourse;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileNameWaterCourse;
  fileControl.ignore(256,'\n');
  ifstream fileWCo(fileNameWaterCourse);
  if (fileWCo == NULL) {
    cout << endl << "Error opening file " << fileNameWaterCourse << endl << endl;
    exit (1);
  }
  // Watercourse elements
  fileWCo.ignore(100,':');
  fileWCo >> numWatc;
  cout << "\n # Number of watercourses " << numWatc << endl;
  WaterCourse *WaterElement = new WaterCourse [numWatc];
  for (i=0; i<numWatc; i++) {
    fileWCo >> j >> ch >> waterCourseId >> lakeNumber >> correction;
    if (j != i) {
      cout << endl << "Error reading file " << fileNameWaterCourse << "\t" << i << "\t" << j << endl;
      exit (1);
    }
    fileWCo.ignore(256,'\n');
    WaterElement[i].SetWaterCourseIndex(i);
    WaterElement[i].SetIdentifier(waterCourseId);
    WaterElement[i].SetLakeNumber(lakeNumber);
    WaterElement[i].SetCorrection(correction);
    WaterElement[i].AllocateAccumulatedDischarge(initialTimeSteps+numberTimeSteps);
    WaterElement[i].AllocateAccumulatedInFlow(initialTimeSteps+numberTimeSteps);
    WaterElement[i].AllocateAccumulatedWaterBalance(initialTimeSteps+numberTimeSteps);
    WaterElement[i].AllocateWaterBalance(initialTimeSteps+numberTimeSteps);
    WaterElement[i].ObsDataInput(startSimulationTime, endSimulationTime, numberTimeSteps, 
                                 ParCommonStore->GetSECONDS_TIMESTEP());
    WaterElement[i].SetSelectedWaterCourseTimeSeriesElements(SelectedWaterCourseTimeSeriesElementsStore);
    cout << " Watercourse index " << i << "  " << "Watercourse identifier " << WaterElement[i].GetIdentifier() << "  " << "Watercourse correction " << WaterElement[i].GetCorrection() << endl;
  }
  // Watercourse outlets
  fileWCo.ignore(100,':');
  fileWCo >> numWatcOut;
  cout << "\n # Number of watercourse outlets " << numWatcOut << endl;
  WaterCourse ** Outlet = new WaterCourse * [numWatcOut];
  for (i=0; i<numWatcOut; i++) {
    fileWCo >> j;
    Outlet[i] = &WaterElement[j];
    cout << " Outlet no. " << i << "\t" << " WaterCourse no. " << j << "\t" << endl;
    fileWCo.ignore(256,'\n');
  }
  // Hierarchy of watercourses
  fileWCo.getline(buffer, 256);
  cout << "\n " << buffer << endl;
  while (fileWCo >> i) {
    fileWCo >> numWatcUp;
    WaterElement[i].SetNumUpStream(numWatcUp);
    fileWCo.ignore(100,':');
    cout << " Downstream, watercourse no.  " << i << "    Identifier  " << WaterElement[i].GetIdentifier() << endl;
    cout << " No. of upstream watercourses " << numWatcUp << endl;
    k = 0;
    while (fileWCo.peek() != '\n') {
      fileWCo >> j;
      cout << "\t" << "Upstream, watercourse no. " << j ;
      WaterElement[i].SetUpStream(k, &WaterElement[j]);
      cout  << "\t" << "UpStream[" << k << "]" << "    Identifier  " << WaterElement[i].GetUpStream(k)->GetIdentifier() << endl;
      while (fileWCo.peek() == ' ') fileWCo.ignore(1,' ');
      k++;
    }
    fileWCo.ignore(256,'\n');
    if (numWatcUp!=k) {
      cout << endl << " Error in number of upstream pointers for watercourse no. " << i << endl << endl;
      exit (1);
    } 
  }
  fileWCo.close();
  // End read file with watercourse information

  // Read information about watercourse elements and landscape elements
  ReadWaterCourseIdentifier(Dew, WaterElement, numWatc, flowHierarchy, fileControl, fout);

  // Read hierarchy of landscape elements
  if (flowHierarchy) {
    ReadLandscapeHierarchy(Dew, fileControl, fout);
  }
  else {
    fileControl.getline(buffer, 256);
    //    fileControl.ignore(256,'\n');
    //    cout << "\n " << buffer << endl;
  }    
  
  // Write information about watercourse outlets and landscape elements to file test_waterland.txt
  //  WriteWaterCourseIdentifier(WaterElement, numWatc, fout);

  // Precipitation and temperature correction for catchments
  int * correctionCatchments = new int[maximumCorrectionCatchments];
  double * correctionPrecipitation = new double[maximumCorrectionCatchments];
  double * correctionTemperature = new double[maximumCorrectionCatchments];
  /*  cout << " File with precipitation and temperature correction for catchments: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  //  cout << fileName << endl;
  ifstream finCorrection(fileName);  // Open for reading
  if (finCorrection == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finCorrection.getline(buffer, 256);
  i=0;
  fout << "Precipitation and temperature correction for catchments: \n";
  while (finCorrection >> correctionCatchments[i]) {
    finCorrection >> correctionPrecipitation[i] >> correctionTemperature[i];
    fout << i << "  " << correctionCatchments[i] << "  " <<  correctionPrecipitation[i] << "  " <<  correctionTemperature[i] << endl;
    i++;
  }
  numberCorrectionCatchments = i;
  for (i=0; i<numWatcOut; i++) { 
    TraverseCorrectionWaterCourse(Outlet[i], numberCorrectionCatchments, correctionCatchments, correctionPrecipitation, correctionTemperature);
  }
  finCorrection.close();
  delete [] correctionCatchments;
  delete [] correctionPrecipitation;
  delete [] correctionTemperature;
  // End precipitation and temperature correction for catchments

  if (inputFormat == 'T' || inputFormat == 't') {
    /*    cout << " File with input data: ";
          cin >> fileName;
          cout << endl;*/
    fileControl.ignore(100,':');
    fileControl >> fileNameInput;
    fileControl.ignore(256,'\n');
    //    cout << fileNameInput << endl;
  }
  else {
    fileControl.getline(buffer, 256);
    //    fileControl.ignore(256,'\n');
    //    cout << "\n " << buffer << endl;
  }    

  // Time series format input data
  if (inputFormat == 'T' || inputFormat == 't') {
    // Object for storing input data for time series
    InputTimeSeries * InputTimeSeriesStore = new
      InputTimeSeries(initialTimeSteps+numberTimeSteps, MetStations->GetNumPrecStations()+MetStations->GetNumTempStations(),
                      startModelTime, endSimulationTime, ParCommonStore->GetSECONDS_TIMESTEP());
    InputTimeSeriesStore->SetCommonPar(ParCommonStore);
    ifstream finInput(fileNameInput);  // Open for reading
    if (finInput == NULL) {
      cout << endl << " Error opening file " << fileNameInput << endl << endl;
      exit(1);
    }
    InputTimeSeriesStore->SetInput(finInput);
    finInput.close();
    InputTimeSeriesStore->WriteInput();

    // Water balance for all elements and time steps for spin-up period
    /*    timeStep=0;
    //    for (time=startModelTime; time<startSimulationTime; time+=1) {
    for (datetime=startModelTime; datetime<startSimulationTime; datetime+=ParCommonStore->GetSECONDS_TIMESTEP()) {
      //    cout << "  " << datetime.getYear() << "  " << datetime.getMonth() << "  " << datetime.getDay() << "  " 
      //         << datetime.getHour() << "  " << datetime.getMinute() << "  " << datetime.getSecond() << endl;

      // ** Algorithm to be performed in case: no input to landscape element from upstream elements
      if (!flowHierarchy) {
        inputDataFound=false;
        firstTotal=true;
        WaterBalanceTimeSeries(Dew, ParCommonStore, MetStations, DateTimeStore,
                               InputTimeSeriesStore, InputElementStore, 
                               numLand, timeStep, initialTimeSteps, &inputDataFound, &firstTotal);
        // Traverse watercourses and landscape elements
        if (inputDataFound) {
          for (i=0; i<numWatcOut; i++) { 
            TraverseWaterCourse(Outlet[i], ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore, 
                                DateTimeStore, timeStep, initialTimeSteps, inputFormat, flowHierarchy, forceDirect, &firstTotal, fout);
          }
        }
        else {
          for (i=0; i<numWatcOut; i++) { 
            TraverseMissingDataWaterCourse(Outlet[i], timeStep, fout);
          }
        }
      }
      // ** End algorithm to be performed in case: no input to landscape element from upstream elements

      // ** Algorithm to be performed in case: input to landscape element from upstream elements
      else {
        firstTotal=true;
        for (i=0; i<numWatcOut; i++) { 
          TraverseWaterCourse(Outlet[i], ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore, 
                              DateTimeStore, timeStep, initialTimeSteps, inputFormat, flowHierarchy, forceDirect, &firstTotal, fout);
        }
      }
      // ** End algorithm to be performed in case: input to landscape element from upstream elements

      timeStep++;
    }
    if (timeStep != initialTimeSteps) {
      cout <<" timeStep != initialTimeSteps " << timeStep << "  " << initialTimeSteps << endl << endl;
      exit(1);
    }*/

    timeStep=0;
    // Water balance for all elements and time steps for spin-up period and simulation period
    //    for (time=startSimulationTime; time<=endSimulationTime; time+=1) {
    //    for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=ParCommonStore->GetSECONDS_TIMESTEP()) {
    for (datetime=startModelTime; datetime<=endSimulationTime; datetime+=ParCommonStore->GetSECONDS_TIMESTEP()) {
      //    cout << "  " << datetime.getYear() << "  " << datetime.getMonth() << "  " << datetime.getDay() << "  " 
      //         << datetime.getHour() << "  " << datetime.getMinute() << "  " << datetime.getSecond() << endl;

      // ** Algorithm to be performed in case: no input to landscape element from upstream elements
      if (!flowHierarchy) {
        inputDataFound=false;
        firstTotal=true;
        WaterBalanceTimeSeries(Dew, ParCommonStore, MetStations, DateTimeStore,
                               InputTimeSeriesStore, InputElementStore, 
                               numLand, timeStep, initialTimeSteps, &inputDataFound, &firstTotal);
        // Traverse watercourses and landscape elements
        if (inputDataFound) {
          for (i=0; i<numWatcOut; i++) { 
            TraverseWaterCourse(Outlet[i], ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore, 
                                DateTimeStore, timeStep, initialTimeSteps, inputFormat, flowHierarchy, forceDirect, &firstTotal, fout);
          }
        }
        else {
          for (i=0; i<numWatcOut; i++) { 
            TraverseMissingDataWaterCourse(Outlet[i], timeStep, fout);
          }
        }
      // Write state variables for all landscape elements 
      /*if (timeStep == (int)((initialTimeSteps+numberTimeSteps-1)/2) || timeStep == initialTimeSteps+numberTimeSteps-1) {
        WriteAsciiGrid(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
        }*/
      }
      // ** End algorithm to be performed in case: no input to landscape element from upstream elements

      // ** Algorithm to be performed in case: input to landscape element from upstream elements
      else {
        firstTotal=true;
        for (i=0; i<numWatcOut; i++) { 
          TraverseWaterCourse(Outlet[i], ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore,
                              DateTimeStore, timeStep, initialTimeSteps, inputFormat, flowHierarchy, forceDirect, &firstTotal, fout);
        }
      }
      // ** End algorithm to be performed in case: input to landscape element from upstream elements

      SnowGlacierIceReDistribution(Outlet, Dew, ParCommonStore, initialTimeSteps, numberTimeSteps, numLand, numWatcOut, timeStep, datetime, 
				   nRows, nCols, noData, xllCorner, yllCorner, cellSize, modelCalibration, fout); 
      timeStep++;
    }
    if (timeStep != initialTimeSteps+numberTimeSteps) {
      cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
      exit(1);
    }
    delete InputTimeSeriesStore;
  }

  // Grid file format input data
  else { 
    // Dummy variable InputTimeSeriesStore
    InputTimeSeries * InputTimeSeriesStore = new
      InputTimeSeries(initialTimeSteps+numberTimeSteps, MetStations->GetNumPrecStations()+MetStations->GetNumTempStations(),
                      startModelTime, endSimulationTime,ParCommonStore->GetSECONDS_TIMESTEP());

    // Path to grid files with meteorological input data
    metPath = getenv("METDATA");
    if (!metPath) {
      cout <<" Environment variable METDATA not defined " << endl << endl;
      exit(1);
    }
    unsigned short int * precip10 = new unsigned short int [nRows*nCols];
    unsigned short int * temp10K = new unsigned short int [nRows*nCols];

    // Water balance for all elements and time steps for spin-up period
    /*    timeStep=0;
    for (datetime=startModelTime; datetime<startSimulationTime; datetime+=ParCommonStore->GetSECONDS_TIMESTEP()) {
      inputDataFound=false;
      firstTotal=true;
      WaterBalanceGrid(Dew, ParCommonStore, InputElementStore, DateTimeStore,
                       numLand, timeStep, initialTimeSteps, nRows, nCols, 
                       datetime, metPath, precip10, temp10K, flowHierarchy, &inputDataFound, &firstTotal);
      // Traverse watercourses and landscape elements
      if (inputDataFound) {
        for (i=0; i<numWatcOut; i++) { 
          //      TraverseWaterCourse(Outlet[i], timeStep, fout);
              TraverseWaterCourse(Outlet[i], ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore, 
                                  DateTimeStore, numberTimeSteps, initialTimeSteps, inputFormat, flowHierarchy, forceDirect, &firstTotal, fout);
        }
      }
      else {
        for (i=0; i<numWatcOut; i++) { 
          TraverseMissingDataWaterCourse(Outlet[i], timeStep, fout);
        }
      }
      timeStep++;
    }
    if (timeStep != initialTimeSteps) {
      cout <<" timeStep != initialTimeSteps " << timeStep << "  " << initialTimeSteps << endl << endl;
      exit(1);
    }*/

    timeStep=0;
    // Water balance for all elements and time steps for spin-up period and simulation period
    //    for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=ParCommonStore->GetSECONDS_TIMESTEP()) {
    for (datetime=startModelTime; datetime<=endSimulationTime; datetime+=ParCommonStore->GetSECONDS_TIMESTEP()) {
      inputDataFound=false;
      firstTotal=true;
      WaterBalanceGrid(Dew, ParCommonStore, InputElementStore, DateTimeStore, 
                       numLand, timeStep, initialTimeSteps, nRows, nCols, 
                       datetime, metPath, precip10, temp10K, flowHierarchy, &inputDataFound, &firstTotal);
      // Traverse watercourses and landscape elements
      if (inputDataFound) {
        for (i=0; i<numWatcOut; i++) { 
          //      TraverseWaterCourse(Outlet[i], timeStep, fout);
              TraverseWaterCourse(Outlet[i], ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore,
                                  DateTimeStore, timeStep, initialTimeSteps, inputFormat, flowHierarchy, forceDirect, &firstTotal, fout);
        }
      }
      else {
        for (i=0; i<numWatcOut; i++) { 
          TraverseMissingDataWaterCourse(Outlet[i], timeStep, fout);
        }
      }
      // Write state variables for all landscape elements 
      /*      if (timeStep == (int)(initialTimeSteps+numberTimeSteps/2.0) || timeStep == initialTimeSteps+numberTimeSteps-1) {*/
      if (timeStep == initialTimeSteps+numberTimeSteps-1) {
        //        WriteReducedBinaryGrid(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
        //        WriteBinaryGrid(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
	//        WriteAsciiGrid(Dew, datetime, numLand, timeStep, nRows, nCols, noData, xllCorner, yllCorner, cellSize, fout);
      }
      SnowGlacierIceReDistribution(Outlet, Dew, ParCommonStore, initialTimeSteps, numberTimeSteps, numLand, numWatcOut, timeStep, datetime, 
				   nRows, nCols, noData, xllCorner, yllCorner, cellSize, modelCalibration, fout); 
      timeStep++;
    }
    if (timeStep != initialTimeSteps+numberTimeSteps) {
      cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
      exit(1);
    }

    delete [] precip10;
    delete [] temp10K;
    delete InputTimeSeriesStore;
  }

  // Write discharge from all watercourse/sub-catchment elements in watercourse hierarchy to output files
  WriteWaterCourseDischarge(WaterElement, numWatc, DateTimeStore, ParCommonStore->GetSECONDS_TIMESTEP(), modelCalibration, fout);
  
  // Write water balance from all watercourse/sub-catchment elements in watercourse hierarchy to output files
  WriteWaterCourseWaterBalance(WaterElement, numWatc, DateTimeStore, ParCommonStore->GetSECONDS_TIMESTEP());

  // Write discharge from all watercourse outlets to output files
  for (i=0; i<numWatcOut; i++)
    WriteOutletDischarge(Outlet[i], startSimulationTime, endSimulationTime, initialTimeSteps, 
                         numberTimeSteps, ParCommonStore->GetSECONDS_TIMESTEP(), false, fout);

  // Write water balance from all watercourse outlets to output files
  for (i=0; i<numWatcOut; i++)
    WriteOutletWaterBalance(Outlet[i], startSimulationTime, endSimulationTime, initialTimeSteps, 
                            numberTimeSteps, ParCommonStore->GetSECONDS_TIMESTEP(), flowHierarchy, forceDirect);

  // Write state variable time series for landscape elements selected for output
  WriteDistributedElementTimeSeries(Dew, numLand, DateTimeStore, ParCommonStore->GetSECONDS_TIMESTEP());  

  // Write total reservoir storage for all time steps
  /*  WriteTotalReservoirStorage(TotalReservoirStore, DateTimeStore, modelCalibration, ParCommonStore->GetSECONDS_TIMESTEP());*/


  /*  k=0;
      fout << "\nPrecipitation correction grid:\n";
      for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {
      if (k<numLand) {
      if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      fout.width(10); 
      fout << Dew[k].GetPrecipitationCorrection() << "  ";
      k++;
      }
      else {
      fout.width(10); fout << noData << "  ";
      }
      }
      else {
      fout.width(10); fout << noData << "  ";
      }
      }
      fout << endl;
      }
      fout << endl;
      k=0;
      fout << "\nTemperature correction grid:\n";
      for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {
      if (k<numLand) {
      if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      fout.width(10); 
      fout << Dew[k].GetTemperatureCorrection() << "  ";
      k++;
      }
      else {
      fout.width(10); fout << noData << "  ";
      }
      }
      else {
      fout.width(10); fout << noData << "  ";
      }
      }
      fout << endl;
      }*/

  
  fout << endl;
  fout.close();
  fileControl.close();

  delete [] WaterElement;
  delete [] Outlet;
  delete [] Dew;
  delete ParCommonStore;
  delete [] ParLandSurfaceStore;
  delete [] ParSubSurfaceHbvStore;
  delete InputElementStore;

  return 0;
}


void SnowGlacierIceReDistribution(WaterCourse ** const Outlet, DistributedElement * const Dew, ParametersCommon * ParCommonStore, int initialTimeSteps, int numberTimeSteps,
                                  int numLand, int numWatcOut, int timeStep, DateTime datetime, 
                                  int nRows, int nCols, int noData, double xllCorner, double yllCorner, double cellSize, bool modelCalibration, ofstream &fout)
{
    int i;
    // Snow store is removed at day no. DAY_SNOW_ZERO
    if (ParCommonStore->GetDAY_SNOW_ZERO() > 0 &&
            dayNumber(datetime.getYear(), datetime.getMonth(), datetime.getDay()) ==
            ParCommonStore->GetDAY_SNOW_ZERO() + leapYear(datetime.getYear()))
    {
      //        cout << "            Snow " << datetime.getYear() << " " << datetime.getMonth() << " " << datetime.getDay() << " " << endl;
        for (i = 0; i < numLand; i++)
        {
            Dew[i].SetSnowStore(0.0);
        }
    }
}


void WaterBalanceTimeSeries(DistributedElement * const Dew, ParametersCommon * const ParCommonStore, 
                            MeteorologicalStations * const MetStations, DateTimeInformation * const DateTimeStore,
                            InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore,
                            int numLand, int timeStep, int initialTimeSteps, bool * inputDataFound, bool * firstTotal)
{
  int i;
  for (i=0; i<numLand; i++) {
    SetDistributedElementInput(&Dew[i], ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore, 
                               timeStep, inputDataFound);
    if (*inputDataFound) {
// ** Algorithm to be performed in case: no input to landscape element from upstream elements
      //      Dew[i].WaterBalance(timeStep);
// ** Algorithm to be performed in case: input to landscape element from upstream elements
      Dew[i].WaterBalance(timeStep,initialTimeSteps,0.0,0.0);
    }
    Dew[i].SetTotalReservoirStorage(DateTimeStore, timeStep, firstTotal, *inputDataFound);
    if (*firstTotal) *firstTotal=false;
  }
}


void SetDistributedElementInput(DistributedElement * const thisElement, ParametersCommon * const ParCommonStore, 
                            MeteorologicalStations * const MetStations, 
                            InputTimeSeries * InputTimeSeriesStore, InputElement * InputElementStore,
                            int timeStep, bool * inputDataFound)
{
  int indSta,indMet,indWeight,indexPrec,indexTemp;
  double precipitation, temperature, weight, elementPrecipitation, elementTemperature;
  double precGradLow, precGradHigh, lapseRate, lapseRateDry, lapseRateWet;
  double elementElevation, gradElevation, metStationElevation, metStationWeight;
  double rainSnowTemperature;
  gradElevation = ParCommonStore->GetGRAD_CHANGE_ALT();
  precGradLow = ParCommonStore->GetPREC_GRAD_LOW();
  precGradHigh = ParCommonStore->GetPREC_GRAD_HIGH();
  lapseRateDry = ParCommonStore->GetLAPSE_DRY();
  lapseRateWet = ParCommonStore->GetLAPSE_WET();
  elementElevation = thisElement->GetElevation();  
  // Precipitation data
  elementPrecipitation = 0.0;
  precipitation = 0.0;
  weight = 0.0;
  for (indexPrec=0; indexPrec<ParCommonStore->GetNUM_PREC_SERIES(); indexPrec++) {
    indSta = thisElement->GetMetSeriesNumber(indexPrec);
    indWeight = indexPrec;
    indMet = indSta;
    metStationElevation = MetStations->GetStationAltitude(indMet);
    metStationWeight = thisElement->GetMetSeriesWeight(indexPrec);
    precipitation = InputTimeSeriesStore->GetInput(timeStep,indMet);
    if (precipitation >= 0.0) { 
      precipitation = precipitation/1000.0;
      if (gradElevation==0.0 || (elementElevation<=gradElevation && metStationElevation<=gradElevation)) 
        precipitation = (precipitation)*pow(precGradLow,(elementElevation-metStationElevation)/100.0)*metStationWeight;
      else if (elementElevation>gradElevation && metStationElevation<gradElevation)
        precipitation = (precipitation)*pow(precGradLow,(gradElevation-metStationElevation)/100.0)*
          pow(precGradHigh,(elementElevation-gradElevation)/100.0)*metStationWeight;
      else if (elementElevation<gradElevation && metStationElevation>gradElevation)
        precipitation = (precipitation)*pow(precGradLow,(elementElevation-gradElevation)/100.0)*
          pow(precGradHigh,(gradElevation-metStationElevation)/100.0)*metStationWeight;
      else
        precipitation = (precipitation)*pow(precGradHigh,(elementElevation-metStationElevation)/100.0)*metStationWeight;
      elementPrecipitation = elementPrecipitation + precipitation;
      weight = weight + metStationWeight;
    }
    /*      cout << indexPrec << "  " << indSta << "  " << indMet << "  " << metStationWeight << "  " << weight << endl;
            if (indexPrec==ParCommonStore->GetNUM_PREC_SERIES()-1) cout << endl;*/
  }
  if (weight > 0.0)
    elementPrecipitation = elementPrecipitation/weight;
  else
    elementPrecipitation = missingData;
  //    cout << indexPrec << "  " << metStationElevation << "  " << weight << "  " << elementPrecipitation << endl;
  /*    for (j=0; j<ParCommonStore->GetNUM_PREC_SERIES(); j++) {
        cout << "P " << thisElement->GetMetSeriesNumber(j) << "  " << thisElement->GetMetSeriesWeight(j) << endl;
        }
        for (j=0; j<ParCommonStore->GetNUM_TEMP_SERIES(); j++) {
        cout << "T " << thisElement->GetMetSeriesNumber(ParCommonStore->GetNUM_PREC_SERIES()+j) << 
        "  " << thisElement->GetMetSeriesWeight(ParCommonStore->GetNUM_PREC_SERIES()+j) << endl;
        }*/
  // Temperature data    
  if (elementPrecipitation == 0.0) 
    lapseRate = lapseRateDry;
  else
    lapseRate = lapseRateWet;
  elementTemperature = 0.0;
  temperature = 0.0;
  weight = 0.0;
  for (indexTemp=0; indexTemp<ParCommonStore->GetNUM_TEMP_SERIES(); indexTemp++) {
    indSta = thisElement->GetMetSeriesNumber(ParCommonStore->GetNUM_PREC_SERIES() + indexTemp);
    indWeight = ParCommonStore->GetNUM_PREC_SERIES() + indexTemp;
    indMet = MetStations->GetNumPrecStations() + indSta;
    metStationElevation = MetStations->GetStationAltitude(indMet);
    metStationWeight = thisElement->GetMetSeriesWeight(indWeight);
    temperature = InputTimeSeriesStore->GetInput(timeStep,indMet); 
    if (temperature > -99.0) { 
      temperature = (temperature/1.0 + lapseRate*(elementElevation-metStationElevation)/100.0)*
        metStationWeight;
      elementTemperature = elementTemperature + temperature; 
      weight = weight + metStationWeight;
    }
    /*      cout << indexTemp << "  " << indSta << "  " << indMet << "  " << metStationWeight << "  " << weight << endl;
            if (indexTemp==ParCommonStore->GetNUM_TEMP_SERIES()-1) cout << endl;*/
  }
  if (weight > 0.0)
    elementTemperature = elementTemperature/weight;
  else
    elementTemperature = missingData;
  if (elementPrecipitation > 0.0) {
    rainSnowTemperature = elementTemperature + 
      lapseRate*(elementElevation-thisElement->GetPrecStationsWeightedElevation())/100.0;
    //  lapseRate*(thisElement->GetTempStationsWeightedElevation()-thisElement->GetPrecStationsWeightedElevation())/100.0;
    if (rainSnowTemperature >= 0.0) 
      elementPrecipitation = elementPrecipitation*ParCommonStore->GetPREC_CORR_RAIN();
    else 
      elementPrecipitation = elementPrecipitation*ParCommonStore->GetPREC_CORR_RAIN()*ParCommonStore->GetPREC_CORR_SNOW();
  }
  //    cout << indexPrec << "  "  << indexTemp << "  " << metStationElevation << "  " << weight << "  " << elementTemperature << endl;
  // Snow store is removed at day no. DAY_SNOW_ZERO
  /*  if (ParCommonStore->GetDAY_SNOW_ZERO() > 0 &&
      dayNumber(InputTimeSeriesStore->GetDateTime(timeStep).getYear(),
                InputTimeSeriesStore->GetDateTime(timeStep).getMonth(),
                InputTimeSeriesStore->GetDateTime(timeStep).getDay()) == 
      ParCommonStore->GetDAY_SNOW_ZERO()+leapYear(InputTimeSeriesStore->GetDateTime(timeStep).getYear())) {
      thisElement->SetSnowStore(0.0);
      //      thisElement->SetSubSurfaceHbvStore(0.2,0.0,0.05);
  }*/
  if (elementPrecipitation > missingData && elementTemperature > missingData) {
    *inputDataFound=true;
    InputElementStore->SetInput(0,elementPrecipitation*thisElement->GetPrecipitationCorrection());
    InputElementStore->SetInput(1,elementTemperature+thisElement->GetTemperatureCorrection());
    //      cout << elementPrecipitation*1000 << "  " << elementTemperature << endl;
    /*      if (dayNumber(InputTimeSeriesStore->GetYear(timeStep),InputTimeSeriesStore->GetMth(timeStep),
            InputTimeSeriesStore->GetDay(timeStep)) == ParCommonStore->GetDAY_SNOW_ZERO()+
            leapYear(InputTimeSeriesStore->GetYear(timeStep))) {
            InputElementStore->SetInput(0,0.0);
            }*/
  }
  else {
    *inputDataFound=false;
    InputElementStore->SetInput(0,missingData);
    InputElementStore->SetInput(1,missingData);
  }
}


void WaterBalanceGrid(DistributedElement * Dew,  ParametersCommon * ParCommonStore, InputElement * InputElementStore,
                      DateTimeInformation * const DateTimeStore,
                      int numLand, int timeStep, int initialTimeSteps, int nRows, int nCols, DateTime datetime, char * metPath,
                      unsigned short int * precip10, unsigned short int * temp10K, bool flowHierarchy, 
                      bool * inputDataFound, bool * firstTotal)
{
  ifstream filePrec, fileTemp;
  int i,j,k;
  unsigned short int metMissing = 10000;
  double precipitation, temperature;
  char fileName[100];
  char precFileName[100];
  char tempFileName[100];
  char hydYear[5];
  strcpy(precFileName,metPath);
  strcpy(tempFileName,metPath);
  strcat(precFileName,"/rr/");
  strcat(tempFileName,"/tm/");
  if (datetime.getMonth() < 9)
    sprintf(hydYear,"%04d",datetime.getYear());
  else
    sprintf(hydYear,"%04d",datetime.getYear()+1);
  strcat(precFileName,hydYear);
  strcat(tempFileName,hydYear);
  sprintf(fileName,"/tm_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());
  strcat(tempFileName,fileName);
  //  if (datetime.getMonth() != 8) {
  sprintf(fileName,"/rr_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());
  /*  }
      else {
      sprintf(fileName,"/rr_%04d_%02d_%02d.bil",datetime.getYear(),7,datetime.getDay());
      }*/
  strcat(precFileName,fileName);
  //  cout << precFileName << "  " << tempFileName << endl;

  filePrec.open(precFileName,ios::in | ios::binary);
  if (filePrec == NULL) {
    cout << endl << "Error opening file " << precFileName << endl << endl;
    exit (1);
  }
  fileTemp.open(tempFileName,ios::in | ios::binary);
  if (fileTemp == NULL) {
    cout << endl << "Error opening file " << tempFileName << endl << endl;
    exit (1);
  }
  //  filePrec.read((unsigned short int*) precip10, sizeof(unsigned short int)*nRows*nCols);
  //  fileTemp.read((unsigned short int*) temp10K, sizeof(unsigned short int)*nRows*nCols);
  streamoff newPosition;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          newPosition = ELEMENT(i,j)*(sizeof(unsigned short int));
          filePrec.seekg(newPosition, ios::beg);
          fileTemp.seekg(newPosition, ios::beg);
          //      filePrec.read((unsigned short int *) &(precip10[ELEMENT(i,j)]), sizeof (unsigned short int));
          //      fileTemp.read((unsigned short int *) &(temp10K[ELEMENT(i,j)]), sizeof (unsigned short int));
          filePrec.read(reinterpret_cast<char *> (&(precip10[ELEMENT(i,j)])), sizeof (unsigned short int));
          fileTemp.read(reinterpret_cast<char *> (&(temp10K[ELEMENT(i,j)])), sizeof (unsigned short int));
          //      printf("%d  %hu  %hu\n",ELEMENT(i,j),precip10[ELEMENT(i,j)],temp10K[ELEMENT(i,j)]);
          if (precip10[ELEMENT(i,j)]<metMissing && temp10K[ELEMENT(i,j)]<metMissing) {
            *inputDataFound=true;
            precipitation = Dew[k].GetPrecipitationCorrection()*(double)precip10[ELEMENT(i,j)]/10000.0;
            temperature = Dew[k].GetTemperatureCorrection()+(double)(temp10K[ELEMENT(i,j)]-2731)/10.0;
            //      printf("%d  %f  %f\n",ELEMENT(i,j),precipitation,temperature);
            InputElementStore->SetInput(0,precipitation);
            InputElementStore->SetInput(1,temperature);
          }
          else {
            //      cout << endl << " Missing meterological data for: " << endl;
            //      cout << precFileName << " or " << tempFileName << endl;
            //      cout << " row = " << i << "  col = " << j << "  element no. = " << ELEMENT(i,j) << endl;
            //      printf("  Precipitation %hu  Temperature %hu\n",precip10,temp10K);
            //      cout << endl << endl;
            //      exit (1);
            //      *inputDataFound=false;
            *inputDataFound=true;
            /*            if (InputElementStore->GetInput(0)==missingData) InputElementStore->SetInput(0,0.0);
                          if (InputElementStore->GetInput(1)==missingData) InputElementStore->SetInput(1,5.0);*/
            InputElementStore->SetInput(0,0.0);
            InputElementStore->SetInput(1,5.0);
          }
          // Snow store is removed at day no. DAY_SNOW_ZERO
	  /*          if (ParCommonStore->GetDAY_SNOW_ZERO() > 0 &&
              dayNumber(datetime.getYear(),datetime.getMonth(),datetime.getDay()) == 
              ParCommonStore->GetDAY_SNOW_ZERO()+leapYear(datetime.getYear())) {
              //      InputElementStore->SetInput(0,0.0);
              Dew[k].SetSnowStore(0.0);
              //      Dew[k].SetSubSurfaceHBVStore(0.2,0.0,0.05);
          }*/
          // Water balace for landscape elements in case: no input to landscape element from upstream elements
          if (!flowHierarchy) {
            Dew[k].WaterBalance(timeStep,initialTimeSteps,0.0,0.0);
            Dew[k].SetTotalReservoirStorage(DateTimeStore, timeStep, firstTotal, *inputDataFound);
            if (*firstTotal) *firstTotal=false;
          }
          // End water balace for landscape elements in case: no input to landscape element from upstream elements
          Dew[k].SetInputValue(0,InputElementStore->GetInput(0));
          Dew[k].SetInputValue(1,InputElementStore->GetInput(1));
          k++;
        }
      }
    }
  }
  filePrec.close();
  fileTemp.close();
}


void ReadWaterCourseIdentifier(DistributedElement * const Dew, WaterCourse * const WaterElement, int numWatc, bool flowHierarchy, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char ch;
  int i,j;
  int waterCourseId, numLandScape;
  int landIndex, geoIndex;
  DistributedElement * thisElement;

  // Read information about watercourse elements and landscape elements
  /*  cout << "\n File with information about watercourse elements and landscape elements: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finWaterCourse(fileName);  // Open for reading
  if (finWaterCourse == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  } 
  // Connect landscape elements to watercourse elements
  for (i=0; i<numWatc; i++) {
    finWaterCourse >> ch >> waterCourseId >> ch >> numLandScape;
    finWaterCourse.ignore(256,'\n');
    WaterElement[i].SetNumLandScape(numLandScape);
    if (waterCourseId!=WaterElement[i].GetIdentifier()) {
      cout << endl << " Error reading file " << fileName << " for watercourse " << i << "\t" 
           << waterCourseId << endl << endl;
      exit (1);
    }
    if (numLandScape > 0) {
      finWaterCourse >> landIndex >> geoIndex;
      thisElement=&Dew[landIndex];
      thisElement->SetLakeNumber(WaterElement[i].GetLakeNumber());
      WaterElement[i].SetLandScapeElement(thisElement);
      //    cout << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
      if (flowHierarchy) {
        thisElement->SetWaterCourseElement(&WaterElement[i]);
      }
      for (j=1; j<numLandScape; j++) {
        finWaterCourse >> landIndex >> geoIndex;
        thisElement->SetNextElement(&Dew[landIndex]);
        thisElement=&Dew[landIndex];
        //      cout << thisElement->GetLandIndex() << "  " << thisElement->GetGeoIndex() << endl;
      }
    }
  }
}


void ReadLandscapeHierarchy(DistributedElement * const Dew, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  int i,j,k;
  int numUp;

  // Read file with pointers between landscape elements
  /*  cout << " File with landscape element hierarchy: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finLandUpFlow(fileName);  // Open for reading
  if (finLandUpFlow == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit (1);
  } 
  while (finLandUpFlow >> i) {
    finLandUpFlow >> numUp;
    Dew[i].SetNumUpLand(numUp);
    finLandUpFlow.ignore(100,':');
    k = 0;
    while (finLandUpFlow.peek() != '\n') {
      finLandUpFlow >> j;
      Dew[i].SetUpLandFlow(k, &Dew[j]);
      k++;
      while (finLandUpFlow.peek() == ' ') finLandUpFlow.ignore(1,' ');
    }
    finLandUpFlow.ignore(256,'\n');
    if (numUp!=k) {
      cout << endl << " Error in number of upland pointers for landscape element no. " 
           << i << endl << endl;
      exit (1);
    } 
  }
  finLandUpFlow.close();
}


// For test purpose only
void WriteWaterCourseIdentifier(WaterCourse * const WaterElement, int numWatc, ofstream &fout)
{
  int i;
  DistributedElement * thisElement;

  // Watercourse outlets
  ofstream waterCourseOut("test_waterland.txt");  // Open for writing"
  for (i=0; i<numWatc; i++) {
    if (WaterElement[i].GetLandScapeElement()) {
      waterCourseOut << "#  " << WaterElement[i].GetIdentifier() << "  #  " 
                   << WaterElement[i].GetNumLandScape() << endl;
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


void TraverseCorrectionWaterCourse(WaterCourse * const thisWaterCourse, int numberCorrectionCatchments,
     int * correctionCatchments, double * correctionPrecipitation, double * correctionTemperature)
{
  int i;
  double precCorr, tempCorr;
  DistributedElement * thisElement;

  for (i=0; i<thisWaterCourse->GetNumUpStream(); i++) {
    TraverseCorrectionWaterCourse(thisWaterCourse->GetUpStream(i), numberCorrectionCatchments,
                                  correctionCatchments, correctionPrecipitation, correctionTemperature);
  }
  if (numberCorrectionCatchments > 0) {
    precCorr = correctionPrecipitation[0];
    tempCorr = correctionTemperature[0];
  }
  else {
    precCorr = 1.0;
    tempCorr = 0.0;
  }
  for (i=0;i<numberCorrectionCatchments;i++) {
    if (thisWaterCourse->GetIdentifier() == correctionCatchments[i]) {
      precCorr = correctionPrecipitation[i];
      tempCorr = correctionTemperature[i];
    }
  }
  thisElement=thisWaterCourse->GetLandScapeElement();
  while (thisElement) {
    TraverseCorrectionLandScape(thisElement, precCorr, tempCorr);
    thisElement = thisElement->GetNextElement();
  }
}


void TraverseCorrectionLandScape(DistributedElement * const thisElement, double precCorr, double tempCorr)
{
  int i;
  for (i=0; i<thisElement->GetNumUpLand(); i++) {
    TraverseCorrectionLandScape(thisElement->GetUpLandFlow(i), precCorr, tempCorr);
  }
  thisElement->SetPrecipitationCorrection(precCorr);
  thisElement->SetTemperatureCorrection(tempCorr);
}


// ** Algorithm to be performed in case: no input to landscape element from upstream elements
//void TraverseWaterCourse(WaterCourse * const thisWaterCourse, int timeStep, ofstream &fout)
// ** Algorithm to be performed in case: input to landscape element from upstream elements
void TraverseWaterCourse(WaterCourse * const thisWaterCourse, ParametersCommon * const ParCommonStore, 
    MeteorologicalStations * const MetStations, InputTimeSeries * InputTimeSeriesStore,
    InputElement * InputElementStore, DateTimeInformation * const DateTimeStore, int timeStep, int initialTimeSteps, char inputFormat,
    bool flowHierarchy, bool forceDirect, bool * firstTotal, ofstream &fout)
{
  int i;
  double accumulatedSum=0.0;
  double accumulatedSumSnow=0.0;
  double accumulatedSumGlacier = 0.0;
  double accumulatedSumHbv=0.0;
  double waterCourseSum=0.0;
  double waterCourseSumSnow=0.0;
  double waterCourseSumGlacier=0.0;
  double waterCourseSumHbv=0.0;
  double accumulatedDischarge=0.0;
  double accumulatedInFlow=0.0;
  double accumulatedPrecipitation=0.0;
  double accumulatedTemperature=0.0;
  double accumulatedSnowStore=0.0;
  double accumulatedGlacierMassBalance = 0.0;
  double accumulatedEvapotranspiration=0.0;
  double accumulatedRunoff=0.0;
  double accumulatedHbvSoilMoisture=0.0;
  double accumulatedHbvSoilMoistureDeficit=0.0;
  double accumulatedHbvPercSoilUpper=0.0;
  double accumulatedHbvUpperZone=0.0;
  double accumulatedHbvLowerZone=0.0;
  double waterCoursePrecipitation=0.0;
  double waterCourseTemperature=0.0;
  double waterCourseSnowStore=0.0;
  double waterCourseGlacierMassBalance = 0.0;
  double waterCourseEvapotranspiration=0.0;
  double waterCourseRunoff=0.0;
  double waterCourseHbvSoilMoisture=0.0;
  double waterCourseHbvSoilMoistureDeficit=0.0;
  double waterCourseHbvPercSoilUpper=0.0;
  double waterCourseHbvUpperZone=0.0;
  double waterCourseHbvLowerZone=0.0;
  bool waterCourseInputFound=true;
  DistributedElement * thisElement;

  //cout << " TraverseWaterCourse " << thisWaterCourse->GetWaterCourseIndex() << "  " << timeStep << endl;
  for (i=0; i<thisWaterCourse->GetNumUpStream(); i++) {
// ** Algorithm to be performed in case: no input to landscape element from upstream elements
    //    TraverseWaterCourse(thisWaterCourse->GetUpStream(i), timeStep, fout);
// ** Algorithm to be performed in case: input to landscape element from upstream elements
    TraverseWaterCourse(thisWaterCourse->GetUpStream(i), ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore, 
                        DateTimeStore, timeStep, initialTimeSteps, inputFormat, flowHierarchy, forceDirect, firstTotal, fout);
    // Fluxes accumulated
    accumulatedDischarge=accumulatedDischarge+thisWaterCourse->GetUpStream(i)->GetAccumulatedDischarge(timeStep);
    accumulatedPrecipitation=accumulatedPrecipitation+thisWaterCourse->GetUpStream(i)->GetAccumulatedPrecipitation(timeStep);
    accumulatedTemperature=accumulatedTemperature+thisWaterCourse->GetUpStream(i)->GetAccumulatedTemperature(timeStep);
    accumulatedEvapotranspiration=accumulatedEvapotranspiration+thisWaterCourse->GetUpStream(i)->GetAccumulatedEvapotranspiration(timeStep);
    accumulatedRunoff=accumulatedRunoff+thisWaterCourse->GetUpStream(i)->GetAccumulatedRunoff(timeStep);
    accumulatedSum=accumulatedSum+thisWaterCourse->GetUpStream(i)->GetAccumulatedSum(timeStep);
    // Snow state variables accumulated
    accumulatedSnowStore=accumulatedSnowStore+thisWaterCourse->GetUpStream(i)->GetAccumulatedSnowStore(timeStep);
    accumulatedSumSnow=accumulatedSumSnow+thisWaterCourse->GetUpStream(i)->GetAccumulatedSumSnow(timeStep);
    // Glacier state variables accumulated
    accumulatedGlacierMassBalance = accumulatedGlacierMassBalance + thisWaterCourse->GetUpStream(i)->GetAccumulatedGlacierMassBalance(timeStep);
    accumulatedSumGlacier = accumulatedSumGlacier + thisWaterCourse->GetUpStream(i)->GetAccumulatedSumGlacier(timeStep);
    // Hbv state variables accumulated
    accumulatedHbvSoilMoisture=accumulatedHbvSoilMoisture+thisWaterCourse->GetUpStream(i)->GetAccumulatedHbvSoilMoisture(timeStep);
    accumulatedHbvSoilMoistureDeficit=accumulatedHbvSoilMoistureDeficit+thisWaterCourse->GetUpStream(i)->GetAccumulatedHbvSoilMoistureDeficit(timeStep);
    accumulatedHbvPercSoilUpper=accumulatedHbvPercSoilUpper+thisWaterCourse->GetUpStream(i)->GetAccumulatedHbvPercSoilUpper(timeStep);
    accumulatedHbvUpperZone=accumulatedHbvUpperZone+thisWaterCourse->GetUpStream(i)->GetAccumulatedHbvUpperZone(timeStep);
    accumulatedHbvLowerZone=accumulatedHbvLowerZone+thisWaterCourse->GetUpStream(i)->GetAccumulatedHbvLowerZone(timeStep);
    accumulatedSumHbv=accumulatedSumHbv+thisWaterCourse->GetUpStream(i)->GetAccumulatedSumHbv(timeStep);
  }
  thisElement=thisWaterCourse->GetLandScapeElement();
  while (thisElement) {
// ** Algorithm to be performed in case: no input to landscape element from upstream elements
    //    TraverseLandScape(thisElement, timeStep, fout);
// ** Algorithm to be performed in case: input to landscape element from upstream elements
    TraverseLandScape(thisElement, ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore, 
                      DateTimeStore, timeStep, initialTimeSteps, inputFormat, &waterCourseInputFound, flowHierarchy, forceDirect, firstTotal, fout);
    // Fluxes accumulated
    accumulatedDischarge=accumulatedDischarge+thisElement->GetAccumulatedDischarge();
    accumulatedInFlow=accumulatedInFlow+thisElement->GetAccumulatedDischarge();
    accumulatedPrecipitation=accumulatedPrecipitation+thisElement->GetAccumulatedPrecipitation();
    accumulatedTemperature=accumulatedTemperature+thisElement->GetAccumulatedTemperature();
    accumulatedEvapotranspiration=accumulatedEvapotranspiration+thisElement->GetAccumulatedEvapotranspiration();
    accumulatedRunoff=accumulatedRunoff+thisElement->GetAccumulatedRunoff();
    accumulatedSum=accumulatedSum+thisElement->GetAccumulatedSum();
    waterCoursePrecipitation=waterCoursePrecipitation+thisElement->GetAccumulatedPrecipitation();
    waterCourseTemperature=waterCourseTemperature+thisElement->GetAccumulatedTemperature();
    waterCourseEvapotranspiration=waterCourseEvapotranspiration+thisElement->GetAccumulatedEvapotranspiration();
    waterCourseRunoff=waterCourseRunoff+thisElement->GetAccumulatedRunoff();
    waterCourseSum=waterCourseSum+thisElement->GetAccumulatedSum();
    // Snow state variables accumulated
    //cout << thisWaterCourse->GetWaterCourseIndex() << "  " << thisElement->GetLandIndex() << "  " << thisElement->GetAccumulatedSumSnow() << "  " << thisElement->GetAccumulatedSnowStore() << endl;
    //cout << waterCourseSumSnow << "  " << waterCourseSnowStore << endl;
    if (thisElement->GetAccumulatedSumSnow() > 0) {
      accumulatedSnowStore=accumulatedSnowStore+thisElement->GetAccumulatedSnowStore();
      accumulatedSumSnow=accumulatedSumSnow+thisElement->GetAccumulatedSumSnow();
      waterCourseSnowStore=waterCourseSnowStore+thisElement->GetAccumulatedSnowStore();
      waterCourseSumSnow=waterCourseSumSnow+thisElement->GetAccumulatedSumSnow();
    }
    //cout << waterCourseSumSnow << "  " << waterCourseSnowStore << endl;
    // Glacier state variables accumulated
    if (thisElement->GetAccumulatedSumGlacier() > 0) {
      accumulatedGlacierMassBalance = accumulatedGlacierMassBalance + thisElement->GetAccumulatedGlacierMassBalance();
      accumulatedSumGlacier = accumulatedSumGlacier + thisElement->GetAccumulatedSumGlacier();
      waterCourseGlacierMassBalance = waterCourseGlacierMassBalance + thisElement->GetAccumulatedGlacierMassBalance();
      waterCourseSumGlacier = waterCourseSumGlacier + thisElement->GetAccumulatedSumGlacier();
    }
    // Hbv state variables accumulated
    if (thisElement->GetAccumulatedSumHbv() > 0) {
      accumulatedHbvSoilMoisture=accumulatedHbvSoilMoisture+thisElement->GetAccumulatedHbvSoilMoisture();
      accumulatedHbvSoilMoistureDeficit=accumulatedHbvSoilMoistureDeficit+thisElement->GetAccumulatedHbvSoilMoistureDeficit();
      accumulatedHbvPercSoilUpper=accumulatedHbvPercSoilUpper+thisElement->GetAccumulatedHbvPercSoilUpper();
      accumulatedHbvUpperZone=accumulatedHbvUpperZone+thisElement->GetAccumulatedHbvUpperZone();
      accumulatedHbvLowerZone=accumulatedHbvLowerZone+thisElement->GetAccumulatedHbvLowerZone();
      accumulatedSumHbv=accumulatedSumHbv+thisElement->GetAccumulatedSumHbv();
      waterCourseHbvSoilMoisture=waterCourseHbvSoilMoisture+thisElement->GetAccumulatedHbvSoilMoisture();
      waterCourseHbvSoilMoistureDeficit=waterCourseHbvSoilMoistureDeficit+thisElement->GetAccumulatedHbvSoilMoistureDeficit();
      waterCourseHbvPercSoilUpper=waterCourseHbvPercSoilUpper+thisElement->GetAccumulatedHbvPercSoilUpper();
      waterCourseHbvUpperZone=waterCourseHbvUpperZone+thisElement->GetAccumulatedHbvUpperZone();
      waterCourseHbvLowerZone=waterCourseHbvLowerZone+thisElement->GetAccumulatedHbvLowerZone();
      waterCourseSumHbv=waterCourseSumHbv+thisElement->GetAccumulatedSumHbv();
    }
    thisElement = thisElement->GetNextElement();
  }
  // ** flowHierarchy == false or true, waterCourseInputFound == true
  if (waterCourseInputFound) {
    // Fluxes accumulated
    thisWaterCourse->SetAccumulatedDischarge(timeStep, accumulatedDischarge);
    thisWaterCourse->SetAccumulatedInFlow(timeStep, accumulatedInFlow);
    if (accumulatedInFlow > accumulatedDischarge) {
      cout << endl << " accumulatedInFlow > accumulatedDischarge  " << accumulatedInFlow << "  " << accumulatedDischarge << endl << endl;
      exit(1);
    }
    thisWaterCourse->SetAccumulatedPrecipitation(timeStep, accumulatedPrecipitation);
    thisWaterCourse->SetAccumulatedTemperature(timeStep, accumulatedTemperature);
    thisWaterCourse->SetAccumulatedEvapotranspiration(timeStep, accumulatedEvapotranspiration);
    thisWaterCourse->SetAccumulatedRunoff(timeStep, accumulatedRunoff);
    thisWaterCourse->SetAccumulatedSum(timeStep, accumulatedSum);
    thisWaterCourse->SetWaterCoursePrecipitation(timeStep, waterCoursePrecipitation/waterCourseSum);
    thisWaterCourse->SetWaterCourseTemperature(timeStep, waterCourseTemperature/waterCourseSum);
    thisWaterCourse->SetWaterCourseEvapotranspiration(timeStep, waterCourseEvapotranspiration/waterCourseSum);
    if (flowHierarchy && !forceDirect) 
      thisWaterCourse->SetWaterCourseRunoff(timeStep, accumulatedDischarge*ParCommonStore->GetSECONDS_TIMESTEP()/waterCourseSum);         /*  discharge (m3/s) -> runoff (m)  */
    else
      thisWaterCourse->SetWaterCourseRunoff(timeStep, waterCourseRunoff/waterCourseSum);
    // Snow state variables accumulated
    //cout << thisWaterCourse->GetWaterCourseIndex() << "  " << waterCourseSumSnow << "  " << waterCourseSnowStore << endl;
    if (waterCourseSumSnow > 0) {
      thisWaterCourse->SetAccumulatedSnowStore(timeStep, accumulatedSnowStore);
      thisWaterCourse->SetAccumulatedSumSnow(timeStep, accumulatedSumSnow);
      thisWaterCourse->SetWaterCourseSnowStore(timeStep, waterCourseSnowStore/waterCourseSumSnow);
    }
    else {
      thisWaterCourse->SetAccumulatedSnowStore(timeStep, missingData);
      thisWaterCourse->SetWaterCourseSnowStore(timeStep, missingData);
    }
    //cout << thisWaterCourse->GetWaterCourseIndex() << "  " << thisWaterCourse->GetWaterCourseSnowStore(timeStep) << endl;
    // Glacier state variables accumulated
    if (waterCourseSumGlacier > 0) {
      thisWaterCourse->SetAccumulatedGlacierMassBalance(timeStep, accumulatedGlacierMassBalance);
      thisWaterCourse->SetAccumulatedSumGlacier(timeStep, accumulatedSumGlacier);
      thisWaterCourse->SetWaterCourseGlacierMassBalance(timeStep, waterCourseGlacierMassBalance/waterCourseSumGlacier);
    }
    else {
      thisWaterCourse->SetAccumulatedGlacierMassBalance(timeStep, accumulatedGlacierMassBalance);
      //      thisWaterCourse->SetAccumulatedSumGlacier(timeStep, accumulatedSumGlacier);
      thisWaterCourse->SetWaterCourseGlacierMassBalance(timeStep, missingData);
    }
    // Hbv state variables accumulated
    if (waterCourseSumHbv > 0) {
      thisWaterCourse->SetAccumulatedHbvSoilMoisture(timeStep, accumulatedHbvSoilMoisture);
      thisWaterCourse->SetAccumulatedHbvSoilMoistureDeficit(timeStep, accumulatedHbvSoilMoistureDeficit);
      thisWaterCourse->SetAccumulatedHbvPercSoilUpper(timeStep, accumulatedHbvPercSoilUpper);
      thisWaterCourse->SetAccumulatedHbvUpperZone(timeStep, accumulatedHbvUpperZone);
      thisWaterCourse->SetAccumulatedHbvLowerZone(timeStep, accumulatedHbvLowerZone);
      thisWaterCourse->SetAccumulatedSumHbv(timeStep, accumulatedSumHbv);
      thisWaterCourse->SetWaterCourseHbvSoilMoisture(timeStep, waterCourseHbvSoilMoisture/waterCourseSumHbv);
      thisWaterCourse->SetWaterCourseHbvSoilMoistureDeficit(timeStep, waterCourseHbvSoilMoistureDeficit/waterCourseSumHbv);
      thisWaterCourse->SetWaterCourseHbvPercSoilUpper(timeStep, waterCourseHbvPercSoilUpper/waterCourseSumHbv);
      thisWaterCourse->SetWaterCourseHbvUpperZone(timeStep, waterCourseHbvUpperZone/waterCourseSumHbv);
      thisWaterCourse->SetWaterCourseHbvLowerZone(timeStep, waterCourseHbvLowerZone/waterCourseSumHbv);
    }
    else {
      thisWaterCourse->SetAccumulatedHbvSoilMoisture(timeStep, missingData);
      thisWaterCourse->SetAccumulatedHbvSoilMoistureDeficit(timeStep, missingData);
      thisWaterCourse->SetAccumulatedHbvPercSoilUpper(timeStep, missingData);
      thisWaterCourse->SetAccumulatedHbvUpperZone(timeStep, missingData);
      thisWaterCourse->SetAccumulatedHbvLowerZone(timeStep, missingData);
      thisWaterCourse->SetWaterCourseHbvSoilMoisture(timeStep, missingData);
      thisWaterCourse->SetWaterCourseHbvSoilMoistureDeficit(timeStep, missingData);
      thisWaterCourse->SetWaterCourseHbvPercSoilUpper(timeStep, missingData);
      thisWaterCourse->SetWaterCourseHbvUpperZone(timeStep, missingData);
      thisWaterCourse->SetWaterCourseHbvLowerZone(timeStep, missingData);
    }
  }
  else {
    SetMissingDataWaterCourse(thisWaterCourse, timeStep, fout);
  }
  //  cout << thisWaterCourse->GetIdentifier() << "  " << timeStep << "  " << "  Discharge: "  
  //       << accumulatedDischarge << endl;
  //       << accumulatedDischarge << "  " << thisWaterCourse->GetAccumulatedDischarge(timeStep) << endl;
}


// ** Algorithm to be performed in case: no input to landscape element from upstream elements
//void TraverseLandScape(DistributedElement * const thisElement, int timeStep, ofstream &fout)
// ** Algorithm to be performed in case: input to landscape element from upstream elements
void TraverseLandScape(DistributedElement * const thisElement, ParametersCommon * const ParCommonStore, 
    MeteorologicalStations * const MetStations, InputTimeSeries * InputTimeSeriesStore,
    InputElement * InputElementStore, DateTimeInformation * const DateTimeStore, int timeStep, int initialTimeSteps, char inputFormat,
    bool * waterCourseInputFound, bool flowHierarchy, bool forceDirect, bool * firstTotal, ofstream &fout)
{
  int i, j;
  double accumulatedSum=0.0;
  double accumulatedSumSnow=0.0;
  double accumulatedSumGlacier = 0.0;
  double accumulatedSumHbv=0.0;
  double accumulatedDischarge=0.0;
  double accumulatedLowerDischarge=0.0;
  double accumulatedUpperDischarge=0.0;
  double accumulatedPrecipitation=0.0;
  double accumulatedTemperature=0.0;
  double accumulatedSnowStore=0.0;
  double accumulatedGlacierMassBalance = 0.0;
  double accumulatedEvapotranspiration=0.0;
  double accumulatedRunoff=0.0;
  double accumulatedHbvSoilMoisture=0.0;
  double accumulatedHbvSoilMoistureDeficit=0.0;
  double accumulatedHbvPercSoilUpper=0.0;
  double accumulatedHbvUpperZone=0.0;
  double accumulatedHbvLowerZone=0.0;
  bool inputDataFound=true;

//  cout << " TraverseLandscape " << thisElement->GetLandIndex() << "  " << timeStep << "  " << thisElement->GetNumUpLand() << endl;
  for (i=0; i<thisElement->GetNumUpLand(); i++) {
// ** Algorithm to be performed in case: no input to landscape element from upstream elements
    //    TraverseLandScape(thisElement->GetUpLandFlow(i), timeStep, fout);
// ** Algorithm to be performed in case: input to landscape element from upstream elements
    TraverseLandScape(thisElement->GetUpLandFlow(i), ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore, 
                      DateTimeStore, timeStep, initialTimeSteps, inputFormat, waterCourseInputFound, flowHierarchy, forceDirect, firstTotal, fout);
    if (thisElement->GetUpLandFlow(i)->GetAccumulatedPrecipitation() != missingData) {
//    if (flowHierarchy && thisElement->GetUpLandFlow(i)->GetAccumulatedPrecipitation() != missingData) {
      // Fluxes accumulated
      accumulatedDischarge=accumulatedDischarge+thisElement->GetUpLandFlow(i)->GetAccumulatedDischarge();
      accumulatedLowerDischarge=accumulatedLowerDischarge+thisElement->GetUpLandFlow(i)->GetAccumulatedLowerDischarge();
      accumulatedUpperDischarge=accumulatedUpperDischarge+thisElement->GetUpLandFlow(i)->GetAccumulatedUpperDischarge();
      accumulatedPrecipitation=accumulatedPrecipitation+thisElement->GetUpLandFlow(i)->GetAccumulatedPrecipitation();
      accumulatedTemperature=accumulatedTemperature+thisElement->GetUpLandFlow(i)->GetAccumulatedTemperature();
      accumulatedEvapotranspiration=accumulatedEvapotranspiration+thisElement->GetUpLandFlow(i)->GetAccumulatedEvapotranspiration();
      accumulatedRunoff=accumulatedRunoff+thisElement->GetUpLandFlow(i)->GetAccumulatedRunoff();
      accumulatedSum=accumulatedSum+thisElement->GetUpLandFlow(i)->GetAccumulatedSum();
      // Snow state variables accumulated
      //cout << thisElement->GetLandIndex() << "  " << accumulatedSnowStore << "  " << thisElement->GetUpLandFlow(i)->GetAccumulatedSnowStore() << "  " << accumulatedSumSnow << endl;
      accumulatedSnowStore=accumulatedSnowStore+thisElement->GetUpLandFlow(i)->GetAccumulatedSnowStore();
      accumulatedSumSnow=accumulatedSumSnow+thisElement->GetUpLandFlow(i)->GetAccumulatedSumSnow();
      //cout << thisElement->GetLandIndex() << "  " << accumulatedSnowStore << "  " << thisElement->GetUpLandFlow(i)->GetAccumulatedSnowStore() << "  " << accumulatedSumSnow << endl;
      // Glacier state variables accumulated
      accumulatedGlacierMassBalance = accumulatedGlacierMassBalance + thisElement->GetUpLandFlow(i)->GetAccumulatedGlacierMassBalance();
      accumulatedSumGlacier = accumulatedSumGlacier + thisElement->GetUpLandFlow(i)->GetAccumulatedSumGlacier();
      // Hbv state variables accumulated
      accumulatedHbvSoilMoisture=accumulatedHbvSoilMoisture+thisElement->GetUpLandFlow(i)->GetAccumulatedHbvSoilMoisture();
      accumulatedHbvSoilMoistureDeficit=accumulatedHbvSoilMoistureDeficit+thisElement->GetUpLandFlow(i)->GetAccumulatedHbvSoilMoistureDeficit();
      accumulatedHbvPercSoilUpper=accumulatedHbvPercSoilUpper+thisElement->GetUpLandFlow(i)->GetAccumulatedHbvPercSoilUpper();
      accumulatedHbvUpperZone=accumulatedHbvUpperZone+thisElement->GetUpLandFlow(i)->GetAccumulatedHbvUpperZone();
      accumulatedHbvLowerZone=accumulatedHbvLowerZone+thisElement->GetUpLandFlow(i)->GetAccumulatedHbvLowerZone();
      accumulatedSumHbv=accumulatedSumHbv+thisElement->GetUpLandFlow(i)->GetAccumulatedSumHbv();
    }
  }

// ** Algorithm to be performed in case: input to landscape element from upstream elements
  if (flowHierarchy) {
    if (inputFormat == 'T' || inputFormat == 't') {
      inputDataFound=false;
      SetDistributedElementInput(thisElement, ParCommonStore, MetStations, InputTimeSeriesStore, InputElementStore, 
                                 timeStep, &inputDataFound);
      if (inputDataFound) {
        if (!forceDirect)
          thisElement->WaterBalance(timeStep, initialTimeSteps, accumulatedLowerDischarge, accumulatedUpperDischarge);
        else
          thisElement->WaterBalance(timeStep, initialTimeSteps, 0.0,0.0);
      }
      else {
        SetMissingDataLandScape(thisElement, timeStep, fout);
        *waterCourseInputFound=false;
      }
      thisElement->SetTotalReservoirStorage(DateTimeStore, timeStep, firstTotal, inputDataFound);
      if (*firstTotal) *firstTotal=false;
    }
    else if (inputFormat == 'G' || inputFormat == 'g') {
      InputElementStore->SetInput(0,thisElement->GetInputValue(0));
      InputElementStore->SetInput(1,thisElement->GetInputValue(1));
      inputDataFound=true;
      if (!forceDirect)
        thisElement->WaterBalance(timeStep,initialTimeSteps,accumulatedLowerDischarge,accumulatedUpperDischarge);
      else
        thisElement->WaterBalance(timeStep,initialTimeSteps,0.0,0.0);
      thisElement->SetTotalReservoirStorage(DateTimeStore, timeStep, firstTotal, inputDataFound);
      if (*firstTotal) *firstTotal=false;
    }
  }
// ** End algorithm to be performed in case: input to landscape element from upstream elements

  // ** flowHierarchy == false or true, inputDataFound == true
  if (inputDataFound) {
    // Fluxes accumulated
    thisElement->SetAccumulatedDischarge(accumulatedDischarge,thisElement->GetDischarge(),flowHierarchy, forceDirect);
    thisElement->SetAccumulatedLowerDischarge(accumulatedLowerDischarge,thisElement->GetLowerDischarge(),flowHierarchy, forceDirect);
    thisElement->SetAccumulatedUpperDischarge(accumulatedUpperDischarge,thisElement->GetUpperDischarge(),flowHierarchy, forceDirect);
    thisElement->SetAccumulatedPrecipitation(accumulatedPrecipitation+thisElement->GetPrecipitation()*thisElement->GetArea());
    thisElement->SetAccumulatedTemperature(accumulatedTemperature+thisElement->GetTemperature()*thisElement->GetArea());
    thisElement->SetAccumulatedEvapotranspiration(accumulatedEvapotranspiration+(thisElement->GetInterceptionLoss()+
                                                  thisElement->GetTranspSoilEvap()+thisElement->GetLakeEvap())*thisElement->GetArea());
    thisElement->SetAccumulatedRunoff(accumulatedRunoff+thisElement->GetRunoff()*thisElement->GetArea());
    thisElement->SetAccumulatedSum(accumulatedSum+thisElement->GetArea());
    //  cout << thisElement->GetAccumulatedSum() << endl;
    // Snow state variables accumulated
    //cout << thisElement->GetLandIndex() << "  " << thisElement->GetSnowStore() << "  " << accumulatedSnowStore << "  " << accumulatedSumSnow << endl;
    if (thisElement->GetSnowStore() != missingData) {
      thisElement->SetAccumulatedSnowStore(accumulatedSnowStore+thisElement->GetSnowStore()*thisElement->GetLandArea());
      thisElement->SetAccumulatedSumSnow(accumulatedSumSnow+thisElement->GetLandArea());
    }
    else {
      thisElement->SetAccumulatedSnowStore(accumulatedSnowStore);
      thisElement->SetAccumulatedSumSnow(accumulatedSumSnow);
        //cout << missingData << "  " << thisElement->GetLandIndex() << "  " << accumulatedSnowStore << "  " << accumulatedSumSnow << endl;
    }
    // Glacier state variables accumulated
    if (thisElement->GetGlacierMassBalance() != missingData)
    {
      thisElement->SetAccumulatedGlacierMassBalance(accumulatedGlacierMassBalance + thisElement->GetGlacierMassBalance()*thisElement->GetArea()*thisElement->GetGlacier()->GetAreaFraction());
      thisElement->SetAccumulatedSumGlacier(accumulatedSumGlacier + thisElement->GetArea()*thisElement->GetGlacier()->GetAreaFraction());
    }
    else {
      thisElement->SetAccumulatedGlacierMassBalance(accumulatedGlacierMassBalance);
      thisElement->SetAccumulatedSumGlacier(accumulatedSumGlacier);
    }
   // Hbv state variables accumulated
    if (thisElement->GetHbvSoilMoisture() != missingData) {
      thisElement->SetAccumulatedHbvSoilMoisture(accumulatedHbvSoilMoisture+thisElement->GetHbvSoilMoisture()*thisElement->GetArea());
      thisElement->SetAccumulatedSumHbv(accumulatedSumHbv+thisElement->GetArea());
    }
    else {
      thisElement->SetAccumulatedHbvSoilMoisture(accumulatedHbvSoilMoisture);
      thisElement->SetAccumulatedSumHbv(accumulatedSumHbv);
    }
    if (thisElement->GetHbvSoilMoistureDeficit() != missingData) {
      thisElement->SetAccumulatedHbvSoilMoistureDeficit(accumulatedHbvSoilMoistureDeficit+thisElement->GetHbvSoilMoistureDeficit()*thisElement->GetArea());
    }
    else {
      thisElement->SetAccumulatedHbvSoilMoistureDeficit(accumulatedHbvSoilMoistureDeficit);
    }
    if (thisElement->GetHbvPercSoilUpper() != missingData) {
      thisElement->SetAccumulatedHbvPercSoilUpper(accumulatedHbvPercSoilUpper+thisElement->GetHbvPercSoilUpper()*thisElement->GetArea());
    }
    else {
      thisElement->SetAccumulatedHbvPercSoilUpper(accumulatedHbvPercSoilUpper);
    }
    if (thisElement->GetHbvUpperZone() != missingData) {
      thisElement->SetAccumulatedHbvUpperZone(accumulatedHbvUpperZone+thisElement->GetHbvUpperZone()*thisElement->GetArea());
    }
    else {
      thisElement->SetAccumulatedHbvUpperZone(accumulatedHbvUpperZone);
    }
    if (thisElement->GetHbvLowerZone() != missingData) {
      thisElement->SetAccumulatedHbvLowerZone(accumulatedHbvLowerZone+thisElement->GetHbvLowerZone()*thisElement->GetArea());
    }
    else {
      thisElement->SetAccumulatedHbvLowerZone(accumulatedHbvLowerZone);
    }
    // Time series for landscape elements
    for (j=0; j<thisElement->GetSelectedHbvTimeSeriesElements()->GetNumberElements(); j++) {
      if (thisElement->GetLandIndex() == thisElement->GetSelectedHbvTimeSeriesElements()->GetHbvTimeSeriesElement(j)) {
        thisElement->SetDistributedElementPrecipitation(timeStep,thisElement->GetPrecipitation());
        thisElement->SetDistributedElementTemperature(timeStep,thisElement->GetTemperature());
        thisElement->SetDistributedElementEvapotranspiration(timeStep,thisElement->GetInterceptionLoss()+
                                                             thisElement->GetTranspSoilEvap()+thisElement->GetLakeEvap());
        thisElement->SetDistributedElementRunoff(timeStep,thisElement->GetRunoff());
        thisElement->SetDistributedElementSnowStore(timeStep,thisElement->GetSnowStore());
        thisElement->SetDistributedElementSnowCoverFraction(timeStep,thisElement->GetSnowCoverFraction());
        thisElement->SetDistributedElementMeltWater(timeStep,thisElement->GetMeltWater());
	thisElement->SetDistributedElementGlacierMassBalance(timeStep, thisElement->GetGlacierMassBalance());
        thisElement->SetDistributedElementHbvSoilMoisture(timeStep,thisElement->GetHbvSoilMoisture());
        thisElement->SetDistributedElementHbvSoilMoistureDeficit(timeStep,thisElement->GetHbvSoilMoistureDeficit());
        thisElement->SetDistributedElementHbvPercSoilUpper(timeStep,thisElement->GetHbvPercSoilUpper());
        thisElement->SetDistributedElementHbvUpperZone(timeStep,thisElement->GetHbvUpperZone());
        thisElement->SetDistributedElementHbvLowerZone(timeStep,thisElement->GetHbvLowerZone());
      }
    }
    for (j=0; j<thisElement->GetSelectedKiWaTimeSeriesElements()->GetNumberElements(); j++) {
      if (thisElement->GetLandIndex() == thisElement->GetSelectedKiWaTimeSeriesElements()->GetKiWaTimeSeriesElement(j)) {
        thisElement->SetDistributedElementPrecipitation(timeStep,thisElement->GetPrecipitation());
        thisElement->SetDistributedElementTemperature(timeStep,thisElement->GetTemperature());
        thisElement->SetDistributedElementEvapotranspiration(timeStep,thisElement->GetInterceptionLoss()+
                     thisElement->GetTranspSoilEvap()+thisElement->GetLakeEvap());
        thisElement->SetDistributedElementRunoff(timeStep,thisElement->GetRunoff());
        thisElement->SetDistributedElementSnowStore(timeStep,thisElement->GetSnowStore());
        thisElement->SetDistributedElementSnowCoverFraction(timeStep,thisElement->GetSnowCoverFraction());
        thisElement->SetDistributedElementMeltWater(timeStep,thisElement->GetMeltWater());
	thisElement->SetDistributedElementGlacierMassBalance(timeStep, thisElement->GetGlacierMassBalance());
        thisElement->SetDistributedElementKiWaSoilMoistureOne(timeStep,thisElement->GetKiWaSoilMoisture
                     (thisElement->GetSelectedKiWaTimeSeriesElements()->GetLengthFractionOne(j)));
        thisElement->SetDistributedElementKiWaSoilMoistureTwo(timeStep,thisElement->GetKiWaSoilMoisture
                     (thisElement->GetSelectedKiWaTimeSeriesElements()->GetLengthFractionTwo(j)));
        thisElement->SetDistributedElementKiWaGroundWaterDepthOne(timeStep,thisElement->GetKiWaGroundWaterDepth
                     (thisElement->GetSelectedKiWaTimeSeriesElements()->GetLengthFractionOne(j)));
        thisElement->SetDistributedElementKiWaGroundWaterDepthTwo(timeStep,thisElement->GetKiWaGroundWaterDepth
                     (thisElement->GetSelectedKiWaTimeSeriesElements()->GetLengthFractionTwo(j)));
      }
    }
  }
}


void TraverseMissingDataWaterCourse(WaterCourse * const thisWaterCourse, int timeStep, ofstream &fout)
{
  int i;
  DistributedElement * thisElement;

  for (i=0; i<thisWaterCourse->GetNumUpStream(); i++) {
    TraverseMissingDataWaterCourse(thisWaterCourse->GetUpStream(i), timeStep, fout);
  }
  thisElement=thisWaterCourse->GetLandScapeElement();
  while (thisElement) {
    SetMissingDataLandScape(thisElement, timeStep, fout);
    thisElement = thisElement->GetNextElement();
  }
  SetMissingDataWaterCourse(thisWaterCourse, timeStep, fout);
}


void SetMissingDataWaterCourse(WaterCourse * const thisWaterCourse, int timeStep, ofstream &fout)
{
  thisWaterCourse->SetAccumulatedDischarge(timeStep, missingData);
  thisWaterCourse->SetAccumulatedInFlow(timeStep, missingData);
  thisWaterCourse->SetAccumulatedPrecipitation(timeStep, missingData);
  thisWaterCourse->SetAccumulatedTemperature(timeStep, missingData);
  thisWaterCourse->SetAccumulatedSnowStore(timeStep, missingData);
  thisWaterCourse->SetAccumulatedGlacierMassBalance(timeStep, missingData);
  thisWaterCourse->SetAccumulatedEvapotranspiration(timeStep, missingData);
  thisWaterCourse->SetAccumulatedRunoff(timeStep, missingData);
  thisWaterCourse->SetAccumulatedHbvSoilMoisture(timeStep, missingData);
  thisWaterCourse->SetAccumulatedHbvSoilMoistureDeficit(timeStep, missingData);
  thisWaterCourse->SetAccumulatedHbvPercSoilUpper(timeStep, missingData);
  thisWaterCourse->SetAccumulatedHbvUpperZone(timeStep, missingData);
  thisWaterCourse->SetAccumulatedHbvLowerZone(timeStep, missingData);
  thisWaterCourse->SetAccumulatedSum(timeStep, missingData);
  thisWaterCourse->SetAccumulatedSumSnow(timeStep, missingData);
  thisWaterCourse->SetAccumulatedSumGlacier(timeStep, missingData);
  thisWaterCourse->SetAccumulatedSumHbv(timeStep, missingData);
  thisWaterCourse->SetWaterCoursePrecipitation(timeStep, missingData);
  thisWaterCourse->SetWaterCourseTemperature(timeStep, missingData);
  thisWaterCourse->SetWaterCourseSnowStore(timeStep, missingData);
  thisWaterCourse->SetWaterCourseGlacierMassBalance(timeStep, missingData);
  thisWaterCourse->SetWaterCourseEvapotranspiration(timeStep, missingData);
  thisWaterCourse->SetWaterCourseRunoff(timeStep, missingData);
  thisWaterCourse->SetWaterCourseHbvSoilMoisture(timeStep, missingData);
  thisWaterCourse->SetWaterCourseHbvSoilMoistureDeficit(timeStep, missingData);
  thisWaterCourse->SetWaterCourseHbvPercSoilUpper(timeStep, missingData);
  thisWaterCourse->SetWaterCourseHbvUpperZone(timeStep, missingData);
  thisWaterCourse->SetWaterCourseHbvLowerZone(timeStep, missingData);
}


void SetMissingDataLandScape(DistributedElement * const thisElement, int timeStep, ofstream &fout)
{
  int j;
  thisElement->SetAccumulatedDischarge(0.0,missingData,false,false);
  thisElement->SetAccumulatedLowerDischarge(0.0,missingData,false,false);
  thisElement->SetAccumulatedUpperDischarge(0.0,missingData,false,false);
  thisElement->SetAccumulatedPrecipitation(missingData);
  thisElement->SetAccumulatedTemperature(missingData);
  thisElement->SetAccumulatedSnowStore(missingData);
  thisElement->SetAccumulatedGlacierMassBalance(missingData);
  thisElement->SetAccumulatedEvapotranspiration(missingData);
  thisElement->SetAccumulatedRunoff(missingData);
  thisElement->SetAccumulatedHbvSoilMoisture(missingData);
  thisElement->SetAccumulatedHbvSoilMoistureDeficit(missingData);
  thisElement->SetAccumulatedHbvPercSoilUpper(missingData);
  thisElement->SetAccumulatedHbvUpperZone(missingData);
  thisElement->SetAccumulatedHbvLowerZone(missingData);
  thisElement->SetAccumulatedSum((int)missingData);
  thisElement->SetAccumulatedSumSnow((int)missingData);
  thisElement->SetAccumulatedSumGlacier(missingData);
  thisElement->SetAccumulatedSumHbv((int)missingData);
  for (j=0; j<thisElement->GetSelectedHbvTimeSeriesElements()->GetNumberElements(); j++) {
    if (thisElement->GetLandIndex() == thisElement->GetSelectedHbvTimeSeriesElements()->GetHbvTimeSeriesElement(j)) {
      thisElement->SetDistributedElementPrecipitation(timeStep,missingData);
      thisElement->SetDistributedElementTemperature(timeStep,missingData);
      thisElement->SetDistributedElementSnowStore(timeStep,missingData);
      thisElement->SetDistributedElementSnowCoverFraction(timeStep,missingData);
      thisElement->SetDistributedElementMeltWater(timeStep,missingData);
      thisElement->SetDistributedElementGlacierMassBalance(timeStep, missingData);
      thisElement->SetDistributedElementEvapotranspiration(timeStep,missingData);
      thisElement->SetDistributedElementRunoff(timeStep,missingData);
      thisElement->SetDistributedElementHbvSoilMoisture(timeStep,missingData);
      thisElement->SetDistributedElementHbvSoilMoistureDeficit(timeStep,missingData);
      thisElement->SetDistributedElementHbvPercSoilUpper(timeStep,missingData);
      thisElement->SetDistributedElementHbvUpperZone(timeStep,missingData);
      thisElement->SetDistributedElementHbvLowerZone(timeStep,missingData);
    }
  }
  for (j=0; j<thisElement->GetSelectedKiWaTimeSeriesElements()->GetNumberElements(); j++) {
    if (thisElement->GetLandIndex() == thisElement->GetSelectedKiWaTimeSeriesElements()->GetKiWaTimeSeriesElement(j)) {
      thisElement->SetDistributedElementPrecipitation(timeStep,missingData);
      thisElement->SetDistributedElementTemperature(timeStep,missingData);
      thisElement->SetDistributedElementSnowStore(timeStep,missingData);
      thisElement->SetDistributedElementSnowCoverFraction(timeStep,missingData);
      thisElement->SetDistributedElementMeltWater(timeStep,missingData);
      thisElement->SetDistributedElementGlacierMassBalance(timeStep, missingData);
      thisElement->SetDistributedElementEvapotranspiration(timeStep,missingData);
      thisElement->SetDistributedElementRunoff(timeStep,missingData);
      thisElement->SetDistributedElementKiWaSoilMoistureOne(timeStep,missingData);
      thisElement->SetDistributedElementKiWaSoilMoistureTwo(timeStep,missingData);
      thisElement->SetDistributedElementKiWaGroundWaterDepthOne(timeStep,missingData);
      thisElement->SetDistributedElementKiWaGroundWaterDepthTwo(timeStep,missingData);
    }
  }
}


//void WriteWaterCourseDischarge(WaterCourse * WaterElement, InputTimeSeries * InputTimeSeriesStore, 
//                              int numWatc, int startSimulationTime, int endSimulationTime, 
//                              int initialTimeSteps, int numberTimeSteps)
void WriteWaterCourseDischarge(WaterCourse * WaterElement, int numWatc, DateTimeInformation * const DateTimeStore,
                                int secondsPerTimeStep, bool modelCalibration, ofstream &fout)
{
  FILE *fpOut, *fpInf;
  char fileName[100];
  int i,j,k,timeStep;
  int initialTimeSteps;
  int numberTimeSteps;
  double sumValue;
  double nashSutcliffe;        /*  Nash-Sutcliffe efficiency criterion  */
  double rootMSError;          /*  Root mean square error criterion  */
  double biasVol;              /*  Bias (volume error) criterion  */
  double pearsCorr;            /*  Pearson's product-moment correlation coefficient  */
  DateTime startModelTime;
  DateTime startSimulationTime;
  DateTime endSimulationTime;
  DateTime datetime;

  initialTimeSteps = DateTimeStore->GetInitialTimeSteps();
  numberTimeSteps = DateTimeStore->GetNumberTimeSteps();
  startModelTime = DateTimeStore->GetStartModelTime();
  startSimulationTime = DateTimeStore->GetStartSimulationTime();
  endSimulationTime = DateTimeStore->GetEndSimulationTime();

  double * observed = new double [numberTimeSteps];
  double * simulated = new double [numberTimeSteps];

  for (i=0; i<numWatc; i++) {
    //    cout << WaterElement[i].GetSelectedWaterCourseTimeSeriesElements()->GetNumberElements() << endl;
    for (k=0; k<WaterElement[i].GetSelectedWaterCourseTimeSeriesElements()->GetNumberElements(); k++) {
      //      cout << WaterElement[i].GetIdentifier() << "    " <<  WaterElement[i].GetSelectedWaterCourseTimeSeriesElements()->GetWaterCourseTimeSeriesElement(k) << endl;
      if (WaterElement[i].GetIdentifier() == WaterElement[i].GetSelectedWaterCourseTimeSeriesElements()->GetWaterCourseTimeSeriesElement(k)) {
        sprintf(fileName,"dew_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpOut = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"inf_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpInf = fopen(fileName, "w")) == NULL ) {
          printf("\n Filen %s ikke funnet!\n\n",fileName);
          exit(1);
        }
        sumValue = 0.0;
        j = 0;
        timeStep=initialTimeSteps;
        for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
          if (WaterElement[i].GetAccumulatedDischarge(timeStep) > missingData) {
            fprintf(fpOut,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetAccumulatedDischarge(timeStep)*WaterElement[i].GetCorrection());
            if (modelCalibration && WaterElement[i].GetObsData(j) > missingData) 
              sumValue = sumValue+WaterElement[i].GetAccumulatedDischarge(timeStep);
          }
          else {
            fprintf(fpOut,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
          }
          if (WaterElement[i].GetAccumulatedDischarge(timeStep) > missingData) {
            fprintf(fpInf,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetAccumulatedInFlow(timeStep));
          }
          else {
            fprintf(fpInf,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
          }
          observed[j] = WaterElement[i].GetObsData(j);
          simulated[j] = WaterElement[i].GetAccumulatedDischarge(timeStep)*WaterElement[i].GetCorrection();
          j++;
          timeStep++;
        }
        if (modelCalibration) fprintf(fpOut,"%29.6f\n", sumValue*WaterElement[i].GetCorrection());
        fclose(fpOut);
        fclose(fpInf);
        if (timeStep != initialTimeSteps+numberTimeSteps) {
          cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
          exit(1);
        }
        ObjectiveCriteria(numberTimeSteps, observed, simulated, &nashSutcliffe, &rootMSError, &biasVol, &pearsCorr);
        cout.precision(2); cout.setf(ios::fixed); cout.setf(ios::showpoint); 
        cout << "\n    Sub-catchment :  " << WaterElement[i].GetIdentifier() << endl << endl;
        cout.width(10); 
        cout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;; 
        cout.width(10); 
        cout << "    Root mean square error                              " << rootMSError << endl;
        cout.width(10); 
        cout << "    Bias (volume error)                                 " << biasVol << endl;
        cout.width(10); 
        cout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
        cout << endl;
        fout.precision(2); fout.setf(ios::fixed); fout.setf(ios::showpoint); 
        fout << "\n    Sub-catchment :  " << WaterElement[i].GetIdentifier() << endl << endl;
        fout.width(10); 
        fout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;; 
        fout.width(10); 
        fout << "    Root mean square error                              " << rootMSError << endl;
        fout.width(10); 
        fout << "    Bias (volume error)                                 " << biasVol << endl;
        fout.width(10); 
        fout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
        fout << endl;
      }
    }
  }

  delete [] observed;
  delete [] simulated;
}


void WriteWaterCourseWaterBalance(WaterCourse * WaterElement, int numWatc, DateTimeInformation * const DateTimeStore,
                                   int secondsPerTimeStep)
{
  FILE *fpPre,*fpTem,*fpSwe, *fpGmb, *fpEva,*fpRun,*fpHsd,*fpHsm,*fpHpe,*fpHuz,*fpHlz,*fpHgw;
  char fileName[100];
  int i,k,timeStep;
  int initialTimeSteps;
  int numberTimeSteps;
  DateTime startModelTime;
  DateTime startSimulationTime;
  DateTime endSimulationTime;
  DateTime datetime;

  initialTimeSteps = DateTimeStore->GetInitialTimeSteps();
  numberTimeSteps = DateTimeStore->GetNumberTimeSteps();
  startModelTime = DateTimeStore->GetStartModelTime();
  startSimulationTime = DateTimeStore->GetStartSimulationTime();
  endSimulationTime = DateTimeStore->GetEndSimulationTime();

  for (i=0; i<numWatc; i++) {
    for (k=0; k<WaterElement[i].GetSelectedWaterCourseTimeSeriesElements()->GetNumberElements(); k++) {
      if (WaterElement[i].GetIdentifier() == WaterElement[i].GetSelectedWaterCourseTimeSeriesElements()->GetWaterCourseTimeSeriesElement(k)) {
        sprintf(fileName,"pre_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpPre = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"tem_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpTem = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"swe_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpSwe = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
	sprintf(fileName, "gmb_%08d.var", WaterElement[i].GetIdentifier());
	if ((fpGmb = fopen(fileName, "w")) == NULL) {
	  printf("\n Filen %s ikke funnet!\n\n", fileName);
	  exit(1);
	}
        sprintf(fileName,"eva_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpEva = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"run_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpRun = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"hsd_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpHsd = fopen(fileName, "w")) == NULL ) {
          printf("\n Filen %s ikke funnet!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"hsm_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpHsm = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"hpe_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpHpe = fopen(fileName, "w")) == NULL ) {
          printf("\n Filen %s ikke funnet!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"huz_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpHuz = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"hlz_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpHlz = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        sprintf(fileName,"hgw_%08d.var",WaterElement[i].GetIdentifier());
        if ((fpHgw = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        timeStep=initialTimeSteps;
        for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
          // Fluxes 
          if (WaterElement[i].GetWaterCoursePrecipitation(timeStep) != missingData) {
            fprintf(fpPre,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCoursePrecipitation(timeStep)*1000.0);
            fprintf(fpTem,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCourseTemperature(timeStep));
            fprintf(fpEva,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCourseEvapotranspiration(timeStep)*1000.0);
            fprintf(fpRun,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCourseRunoff(timeStep)*1000.0);
          }
          else {
            fprintf(fpPre,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
            fprintf(fpTem,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
            fprintf(fpEva,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
            fprintf(fpRun,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
          }
          // Snow state variables
          if (WaterElement[i].GetWaterCourseSnowStore(timeStep) != missingData) 
            fprintf(fpSwe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCourseSnowStore(timeStep)*1000.0);
          else
            fprintf(fpSwe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
 
	  // Glacier amass balance
	  if (WaterElement[i].GetWaterCourseGlacierMassBalance(timeStep) != missingData) {
	    fprintf(fpGmb, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
		    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
		    WaterElement[i].GetWaterCourseGlacierMassBalance(timeStep) * 1000.0);
          }
          else {
	    fprintf(fpGmb, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
		    datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
	  }
	  // Hbv state variables
          if (WaterElement[i].GetWaterCourseHbvSoilMoisture(timeStep) != missingData) { 
            fprintf(fpHsd,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCourseHbvSoilMoistureDeficit(timeStep)*1000.0);
            fprintf(fpHsm,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCourseHbvSoilMoisture(timeStep)*1000.0);
            fprintf(fpHpe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCourseHbvPercSoilUpper(timeStep)*1000.0);
            fprintf(fpHuz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCourseHbvUpperZone(timeStep)*1000.0);
            fprintf(fpHlz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    WaterElement[i].GetWaterCourseHbvLowerZone(timeStep)*1000.0);
            fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    (WaterElement[i].GetWaterCourseHbvUpperZone(timeStep)+WaterElement[i].GetWaterCourseHbvLowerZone(timeStep))*1000.0);
          }
          else {
            fprintf(fpHsd,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
            fprintf(fpHsm,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
            fprintf(fpHpe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
            fprintf(fpHuz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
            fprintf(fpHlz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
            fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
          }
          timeStep++;
        }
        fclose(fpPre);
        fclose(fpTem);
        fclose(fpSwe);
        fclose(fpEva);
        fclose(fpRun);
        fclose(fpHsd);
        fclose(fpHsm);
        fclose(fpHpe);
        fclose(fpHuz);
        fclose(fpHlz);
        fclose(fpHgw);
        if (timeStep != initialTimeSteps+numberTimeSteps) {
          cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
          exit(1);
        }
      }
    }
  }
}


void WriteOutletDischarge(WaterCourse * const WaterElement, DateTime startSimulationTime, 
                          DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                          int secondsPerTimeStep, bool modelCalibration, ofstream &fout)
{
  FILE *fpOut;
  char fileName[100];
  double sumValue;
  double nashSutcliffe;        /*  Nash-Sutcliffe efficiency criterion  */
  double rootMSError;          /*  Root mean square error criterion  */
  double biasVol;              /*  Bias (volume error) criterion  */
  double pearsCorr;            /*  Pearson's product-moment correlation coefficient  */
  int i,j,k,timeStep;
  DateTime datetime;
  double * observed = new double [numberTimeSteps];
  double * simulated = new double [numberTimeSteps];

  sprintf(fileName,"dew_%08d.out",WaterElement->GetIdentifier());
  if ((fpOut = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sumValue = 0.0;
  j = 0;
  timeStep=initialTimeSteps;
  for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
    if (WaterElement->GetAccumulatedDischarge(timeStep) > missingData) {
      fprintf(fpOut,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedDischarge(timeStep)*WaterElement->GetCorrection());
      if (modelCalibration && WaterElement->GetObsData(j) > missingData) 
        sumValue = sumValue+WaterElement->GetAccumulatedDischarge(timeStep);
    }
    else {
      fprintf(fpOut,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
              datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
    }
    observed[j] = WaterElement->GetObsData(j);
    simulated[j] = WaterElement->GetAccumulatedDischarge(timeStep)*WaterElement->GetCorrection();
    j++;
    timeStep++;
  }
  if (modelCalibration) fprintf(fpOut,"%29.6f\n", sumValue*WaterElement->GetCorrection());
  fclose(fpOut);
  if (timeStep != initialTimeSteps+numberTimeSteps) {
    cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
    exit(1);
  }
  ObjectiveCriteria(numberTimeSteps, observed, simulated, &nashSutcliffe, &rootMSError, &biasVol, &pearsCorr);
  cout.precision(2); cout.setf(ios::fixed); cout.setf(ios::showpoint); 
  cout << "\n    Watercourse outlet :  " << WaterElement->GetIdentifier() << endl << endl;
  cout.width(10); 
  cout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;; 
  cout.width(10); 
  cout << "    Root mean square error                              " << rootMSError << endl;
  cout.width(10); 
  cout << "    Bias (volume error)                                 " << biasVol << endl;
  cout.width(10); 
  cout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
  cout << endl;
  fout.precision(2); fout.setf(ios::fixed); fout.setf(ios::showpoint); 
  fout << "\n    Watercourse outlet :  " << WaterElement->GetIdentifier() << endl << endl;
  fout.width(10); 
  fout << "    Nash-Sutcliffe efficiency                           " << nashSutcliffe << endl;; 
  fout.width(10); 
  fout << "    Root mean square error                              " << rootMSError << endl;
  fout.width(10); 
  fout << "    Bias (volume error)                                 " << biasVol << endl;
  fout.width(10); 
  fout << "    Pearson's product-moment correlation coefficient    " << pearsCorr << endl;
  fout << endl;

  delete [] observed;
  delete [] simulated;
}


void WriteOutletWaterBalance(WaterCourse * const WaterElement, DateTime startSimulationTime, 
                             DateTime endSimulationTime, int initialTimeSteps, int numberTimeSteps,
                             int secondsPerTimeStep, bool flowHierarchy, bool forceDirect)
{
  FILE *fpPre,*fpTem,*fpSwe,*fpGmb,*fpEva,*fpRun,*fpHsd,*fpHsm,*fpHpe,*fpHuz,*fpHlz,*fpHgw;
  char fileName[100];
  int i,k,timeStep;
  DateTime datetime;

  sprintf(fileName,"pre_%08d.out",WaterElement->GetIdentifier());
  if ((fpPre = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName,"tem_%08d.out",WaterElement->GetIdentifier());
  if ((fpTem = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName,"swe_%08d.out",WaterElement->GetIdentifier());
  if ((fpSwe = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName, "gmb_%08d.out", WaterElement->GetIdentifier());
  if ((fpGmb = fopen(fileName, "w")) == NULL) {
    printf("\n Filen %s ikke funnet!\n\n", fileName);
    exit(1);
  }
  sprintf(fileName,"eva_%08d.out",WaterElement->GetIdentifier());
  if ((fpEva = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName,"run_%08d.out",WaterElement->GetIdentifier());
  if ((fpRun = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName,"hsd_%08d.out",WaterElement->GetIdentifier());
  if ((fpHsd = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName,"hsm_%08d.out",WaterElement->GetIdentifier());
  if ((fpHsm = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName,"hpe_%08d.out",WaterElement->GetIdentifier());
  if ((fpHpe = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName,"huz_%08d.out",WaterElement->GetIdentifier());
  if ((fpHuz = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName,"hlz_%08d.out",WaterElement->GetIdentifier());
  if ((fpHlz = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }
  sprintf(fileName,"hgw_%08d.out",WaterElement->GetIdentifier());
  if ((fpHgw = fopen(fileName, "w")) == NULL ) {
    printf("\n Filen %s ikke funnet!\n\n",fileName);
    exit(1);
  }

  timeStep=initialTimeSteps;
  for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
    // Fluxes 
    if (WaterElement->GetAccumulatedPrecipitation(timeStep) != missingData) {
      fprintf(fpPre,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedPrecipitation(timeStep)*1000.0/WaterElement->GetAccumulatedSum(timeStep));
      fprintf(fpTem,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedTemperature(timeStep)/WaterElement->GetAccumulatedSum(timeStep));
      fprintf(fpEva,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedEvapotranspiration(timeStep)*1000.0/WaterElement->GetAccumulatedSum(timeStep));
      if (flowHierarchy && !forceDirect) {
        fprintf(fpRun,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                WaterElement->GetAccumulatedDischarge(timeStep)*secondsPerTimeStep*1000.0/WaterElement->GetAccumulatedSum(timeStep));
      }
      else {
        fprintf(fpRun,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                WaterElement->GetAccumulatedRunoff(timeStep)*1000.0/WaterElement->GetAccumulatedSum(timeStep));
      }
    }
    else {
      fprintf(fpPre,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),  
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      fprintf(fpTem,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      fprintf(fpEva,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      fprintf(fpRun,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
    }
    // Snow state variables
    if (WaterElement->GetAccumulatedSnowStore(timeStep) != missingData) 
      fprintf(fpSwe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedSnowStore(timeStep)*1000.0/WaterElement->GetAccumulatedSumSnow(timeStep));
    else
      fprintf(fpSwe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
    // Glacier mass balance
    if (WaterElement->GetAccumulatedGlacierMassBalance(timeStep) != missingData && WaterElement->GetAccumulatedSumGlacier(timeStep) > 0.0) {
      fprintf(fpGmb, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
	      datetime.getDay(), datetime.getHour(), datetime.getMinute(),
	      WaterElement->GetAccumulatedGlacierMassBalance(timeStep) * 1000.0 / WaterElement->GetAccumulatedSumGlacier(timeStep));
    }
    else {
      fprintf(fpGmb, "%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(),
	      datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
    }
    // Hbv state variables
    if (WaterElement->GetAccumulatedHbvSoilMoisture(timeStep) != missingData) { 
      fprintf(fpHsd,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedHbvSoilMoistureDeficit(timeStep)*1000.0/WaterElement->GetAccumulatedSumHbv(timeStep));
      fprintf(fpHsm,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedHbvSoilMoisture(timeStep)*1000.0/WaterElement->GetAccumulatedSumHbv(timeStep));
      fprintf(fpHpe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedHbvPercSoilUpper(timeStep)*1000.0/WaterElement->GetAccumulatedSumHbv(timeStep));
      fprintf(fpHuz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedHbvUpperZone(timeStep)*1000.0/WaterElement->GetAccumulatedSumHbv(timeStep));
      fprintf(fpHlz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              WaterElement->GetAccumulatedHbvLowerZone(timeStep)*1000.0/WaterElement->GetAccumulatedSumHbv(timeStep));
      fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              (WaterElement->GetAccumulatedHbvUpperZone(timeStep)+WaterElement->GetAccumulatedHbvLowerZone(timeStep))*1000.0/WaterElement->GetAccumulatedSumHbv(timeStep));
    }
    else {
      fprintf(fpHsd,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      fprintf(fpHsm,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      fprintf(fpHpe,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      fprintf(fpHuz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      fprintf(fpHlz,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
      fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
    }
    timeStep++;
  }
  fclose(fpPre);
  fclose(fpTem);
  fclose(fpSwe);
  fclose(fpEva);
  fclose(fpRun);
  fclose(fpHsd);
  fclose(fpHsm);
  fclose(fpHpe);
  fclose(fpHuz);
  fclose(fpHlz);
  if (timeStep != initialTimeSteps+numberTimeSteps) {
    cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
    exit(1);
  }
}


void WriteDistributedElementTimeSeries(DistributedElement * const Dew, int numLand, DateTimeInformation * const DateTimeStore,
                                       int secondsPerTimeStep)
{
  FILE *fpHgw, *fpKgw1, *fpKgw2;
  char fileName[100];
  int i,j,timeStep;
  int initialTimeSteps;
  int numberTimeSteps;
  DateTime startModelTime;
  DateTime startSimulationTime;
  DateTime endSimulationTime;
  DateTime datetime;

  initialTimeSteps = DateTimeStore->GetInitialTimeSteps();
  numberTimeSteps = DateTimeStore->GetNumberTimeSteps();
  startModelTime = DateTimeStore->GetStartModelTime();
  startSimulationTime = DateTimeStore->GetStartSimulationTime();
  endSimulationTime = DateTimeStore->GetEndSimulationTime();

  for (i=0; i<numLand; i++) {
    //  HBV state variables time series output
    for (j=0; j<Dew[i].GetSelectedHbvTimeSeriesElements()->GetNumberElements(); j++) {
      if (Dew[i].GetLandIndex() == Dew[i].GetSelectedHbvTimeSeriesElements()->GetHbvTimeSeriesElement(j)) {
        sprintf(fileName,"HBV_groundwater_%d.var",Dew[i].GetLandIndex());
        if ((fpHgw = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        timeStep=initialTimeSteps;
        for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
          if (Dew[i].GetDistributedElementHbvUpperZone(timeStep) != missingData) { 
            fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    Dew[i].GetSelectedHbvTimeSeriesElements()->GetGroundWaterRef(j) + 
                    (Dew[i].GetDistributedElementHbvUpperZone(timeStep)+Dew[i].GetDistributedElementHbvLowerZone(timeStep)) /
                    Dew[i].GetSelectedHbvTimeSeriesElements()->GetEffPor(j));
          }
          else {
            fprintf(fpHgw,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
          }
          timeStep++;
        }
        fclose(fpHgw);
        if (timeStep != initialTimeSteps+numberTimeSteps) {
          cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
          exit(1);
        }
      }
    }
    //  KiWa state variables time series output
    for (j=0; j<Dew[i].GetSelectedKiWaTimeSeriesElements()->GetNumberElements(); j++) {
      if (Dew[i].GetLandIndex() == Dew[i].GetSelectedKiWaTimeSeriesElements()->GetKiWaTimeSeriesElement(j)) {
        sprintf(fileName,"KiWa_groundwater_%d_one.var",Dew[i].GetLandIndex());
        if ((fpKgw1 = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        timeStep=initialTimeSteps;
        for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
          if (Dew[i].GetDistributedElementKiWaGroundWaterDepthOne(timeStep) != missingData) { 
            fprintf(fpKgw1,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    Dew[i].GetDistributedElementKiWaGroundWaterDepthOne(timeStep));
          }
          else {
            fprintf(fpKgw1,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
          }
          timeStep++;
        }
        fclose(fpKgw1);
        if (timeStep != initialTimeSteps+numberTimeSteps) {
          cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
          exit(1);
        }
        sprintf(fileName,"KiWa_groundwater_%d_two.var",Dew[i].GetLandIndex());
        if ((fpKgw2 = fopen(fileName, "w")) == NULL ) {
          printf("\n File %s not found!\n\n",fileName);
          exit(1);
        }
        timeStep=initialTimeSteps;
        for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
          if (Dew[i].GetDistributedElementKiWaGroundWaterDepthTwo(timeStep) != missingData) { 
            fprintf(fpKgw2,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),
                    Dew[i].GetDistributedElementKiWaGroundWaterDepthTwo(timeStep));
          }
          else {
            fprintf(fpKgw2,"%4d%02d%02d/%02d%02d%16.6f\n", datetime.getYear(), datetime.getMonth(), 
                    datetime.getDay(), datetime.getHour(), datetime.getMinute(),missingData);
          }
          timeStep++;
        }
        fclose(fpKgw2);
        if (timeStep != initialTimeSteps+numberTimeSteps) {
          cout <<" timeStep != initialTimeSteps+numberTimeSteps " << timeStep << "  " << initialTimeSteps+numberTimeSteps << endl << endl;
          exit(1);
        }
      }
    }
  }
}


void WriteTotalReservoirStorage(TotalReservoirStorage * const TotalReservoirStore, DateTimeInformation * const DateTimeStore,
                                bool modelCalibration, int secondsPerTimeStep)
{
  FILE *fpOut;
  char fileName[100];
  int timeStep;
  int initialTimeSteps;
  int numberTimeSteps;
  DateTime startModelTime;
  DateTime startSimulationTime;
  DateTime endSimulationTime;
  DateTime datetime;

  initialTimeSteps = DateTimeStore->GetInitialTimeSteps();
  numberTimeSteps = DateTimeStore->GetNumberTimeSteps();
  startModelTime = DateTimeStore->GetStartModelTime();
  startSimulationTime = DateTimeStore->GetStartSimulationTime();
  endSimulationTime = DateTimeStore->GetEndSimulationTime();

  sprintf(fileName,"res_totalvol.var");
  if ((fpOut = fopen(fileName, "w")) == NULL ) {
      printf("\n File %s not found!\n\n",fileName);
      exit(1);
  }
  if (modelCalibration) 
      fprintf(fpOut,"%29.6f\n", TotalReservoirStore->GetTotalReservoirStorage(initialTimeSteps+numberTimeSteps-1) - 
              TotalReservoirStore->GetTotalReservoirStorage(initialTimeSteps));
  timeStep=initialTimeSteps;
  for (datetime=startSimulationTime; datetime<=endSimulationTime; datetime+=secondsPerTimeStep) {
    if (TotalReservoirStore->GetTotalReservoirStorage(timeStep) != missingData) {
      fprintf(fpOut,"%4d%02d%02d/%02d%02d%29.6f\n", datetime.getYear(), datetime.getMonth(),  
              datetime.getDay(), datetime.getHour(), datetime.getMinute(),
              TotalReservoirStore->GetTotalReservoirStorage(timeStep));
    }
    else {
      fprintf(fpOut,"%4d%02d%02d/%02d%02d%29.6f\n", datetime.getYear(), datetime.getMonth(),  
              datetime.getDay(), datetime.getHour(), datetime.getMinute(), missingData);
    }
    timeStep++;
  }
  fclose(fpOut);
}


/*  Objective criteria of model performance  */
void ObjectiveCriteria(int numberTimeSteps, double *obs, double *sim, 
                        double *ns, double *rmse, double *bias, double *pears)
{
  int i, sum_count;
  double obs_mean, sim_mean;
  double sum_1;
  double sum_2;
  double sum_3;

  *ns = missingData;            /*  Nash-Sutcliffe efficiency criterion  */
  *rmse = missingData;          /*  Root mean square error criterion  */
  *bias = missingData;          /*  Bias (volume error) criterion  */
  *pears = missingData;         /*  Pearson's product-moment correlation coefficient  */
  /*    for (i = 0; i < numberTimeSteps; i++) sim[i] = obs[i];*/

  /* Mean values */
  obs_mean = 0; 
  sim_mean = 0;
  sum_count = 0;
  for (i = 0; i < numberTimeSteps; i++) {
    if (sim[i] > missingData && obs[i] > missingData) {
      obs_mean = obs_mean + obs[i];
      sim_mean = sim_mean + sim[i];
      sum_count++;
    }
  }
  if (sum_count > 0) {
    obs_mean = obs_mean / (double)sum_count;
    sim_mean = sim_mean / (double)sum_count;
    /* Nash-Sutcliffe efficiency criterion */
    sum_1 = 0;
    sum_2 = 0;
    for (i = 0; i < numberTimeSteps; i++) {
      if (sim[i] > missingData && obs[i] > missingData) {
        sum_1 = sum_1 + pow((obs[i] - sim[i]),2);
        sum_2 = sum_2 + pow((obs[i] - obs_mean),2);
      }
    }
    *ns = 1.0-(sum_1/sum_2);
    /* Root mean square error criterion */
    sum_1 = 0;
    for (i = 0; i < numberTimeSteps; i++) {
      if (sim[i] > missingData && obs[i] > missingData) {
        sum_1 = sum_1 + pow((obs[i] - sim[i]),2);
      }
    }
    *rmse = pow((sum_1/(double)i),0.5);
    /* Bias (volume error) criterion */
    sum_1 = 0;
    sum_2 = 0;
    for (i = 0; i < numberTimeSteps; i++) {
      if (sim[i] > missingData && obs[i] > missingData) {
        sum_1 = sum_1 + (sim[i] - obs[i]);
        sum_2 = sum_2 + obs[i];
      }
    }
    *bias = sum_1/sum_2;
    /* Pearson's product-moment correlation coefficient criterion */
    sum_1 = 0;
    sum_2 = 0;
    sum_3 = 0;
    for (i = 0; i < numberTimeSteps; i++) {
      if (sim[i] > missingData && obs[i] > missingData) {
        sum_1 = sum_1 + (obs[i] - obs_mean) * (sim[i] - sim_mean);
        sum_2 = sum_2 + pow((obs[i] - obs_mean),2);
        sum_3 = sum_3 + pow((sim[i] - sim_mean),2);
      }
    }
    *pears = pow((sum_1 / sqrt(sum_2 * sum_3)),2);
  }    
  return;
}


void WriteReducedBinaryGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,  
                     int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
  int i,j,k;
  char fileName[100];
  char dirName[100];
  char timeName[100];
  //  unsigned short int noData16=65535;
  ofstream fileTemp, filePre, fileEva, fileSwe, fileScf, fileSmw, fileRun, fileDew, fileHsd, fileHsm, fileHgw, fileHlz;
  sprintf(dirName,"./%04d/",datetime.getYear());
  sprintf(timeName,"_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());

  //  Temperature in landscape elements 
  float * temp = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"tem");
  strcat(fileName,timeName);
  fileTemp.open(fileName,ios::out | ios::binary);
  if (fileTemp == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//    if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
//      if (Dew[k].GetAccumulatedTemperature() > missingData) 
      if (Dew[k].GetTemperature() > missingData) 
        temp[k] = (float) ((Dew[k].GetTemperature()));
      else
        temp[k] = (float) noData;
      k++;
  }
/*        else {
          temp[ELEMENT(i,j)] = noData;
        }
      }
      else {
        temp[ELEMENT(i,j)] = noData;
      }
      }
      }*/
  fileTemp.write(reinterpret_cast<char *>(temp), sizeof(float)*numLand);
  fileTemp.close();
  delete [] temp;
  
  //  Precipitation in landscape elements 
  float * pre = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"pre");
  strcat(fileName,timeName);
  filePre.open(fileName,ios::out | ios::binary);
  if (filePre == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//    if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
//      if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
      if (Dew[k].GetPrecipitation() > missingData) 
        pre[k] = (float) ((Dew[k].GetPrecipitation())*1000.0);
      else
        pre[k] = (float) noData;
      k++;
  }
/*        else {
          pre[ELEMENT(i,j)] = noData;
          }
          }
          else {
          pre[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  filePre.write(reinterpret_cast<char *>(pre), sizeof(float)*numLand);
  filePre.close();
  delete [] pre;
  
  //  Evapotranspiration from landscape elements
  float * eva = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"eva");
  strcat(fileName,timeName);
  fileEva.open(fileName,ios::out | ios::binary);
  if (fileEva == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//      if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
    if (Dew[k].GetInterceptionLoss() > missingData &&
        Dew[k].GetTranspSoilEvap() > missingData &&
        Dew[k].GetLakeEvap() > missingData)
      eva[k] = (float) ((Dew[k].GetInterceptionLoss() + Dew[k].GetTranspSoilEvap() + 
                         Dew[k].GetLakeEvap())*1000.0);
    else
      eva[k] = (float) noData;
    k++;
  }
/*        else {
          eva[ELEMENT(i,j)] = noData;
          }
          }
          else {
          eva[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileEva.write(reinterpret_cast<char *>(eva), sizeof(float)*numLand);
  fileEva.close();
  delete [] eva;
  
  //  Snow store in landscape elements
  float * swe = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"swe");
  strcat(fileName,timeName);
  fileSwe.open(fileName,ios::out | ios::binary);
  if (fileSwe == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (Dew[k].GetSnowStore() > missingData) 
          swe[k] = (float) ((Dew[k].GetSnowStore())*1000.0);
      else
          swe[k] = (float) noData;
      k++;
  }
/*        else {
          swe[ELEMENT(i,j)] = noData;
          }
          }
          else {
          swe[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileSwe.write(reinterpret_cast<char *>(swe), sizeof(float)*numLand);
  fileSwe.close();
  delete [] swe;
  
  //  Snowcover fraction in landscape elements
  float * scf = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"scf");
  strcat(fileName,timeName);
  fileScf.open(fileName,ios::out | ios::binary);
  if (fileScf == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (Dew[k].GetSnowCoverFraction() > missingData) 
          scf[k] = (float) ((Dew[k].GetSnowCoverFraction())*100.0);
      else
          scf[k] = (float) noData;
      k++;
  }
/*        else {
          scf[ELEMENT(i,j)] = noData;
          }
          }
          else {
          scf[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileScf.write(reinterpret_cast<char *>(scf), sizeof(float)*numLand);
  fileScf.close();
  delete [] scf;
   
  //  Snow meltwater in landscape elements
  float * smw = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"smw");
  strcat(fileName,timeName);
  fileSmw.open(fileName,ios::out | ios::binary);
  if (fileSmw == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (Dew[k].GetMeltWater() > missingData) 
          smw[k] = (float) ((Dew[k].GetMeltWater())*1000.0);
      else
          smw[k] = (float) noData;
      k++;
  }
/*        else {
          smw[ELEMENT(i,j)] = noData;
          }
          }
          else {
          smw[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileSmw.write(reinterpret_cast<char *>(smw), sizeof(float)*numLand);
  fileSmw.close();
  delete [] smw;
 
  //  Runoff from landscape elements
  float * runoff = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"run");
  strcat(fileName,timeName);
  fileRun.open(fileName,ios::out | ios::binary);
  if (fileRun == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (Dew[k].GetRunoff() > missingData) 
          runoff[k] = (float) ((Dew[k].GetRunoff())*1000.0);
      else
          runoff[k] = (float) noData;
      k++;
  }
/*        else {
          runoff[ELEMENT(i,j)] = noData;
          }
          }
          else {
          runoff[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileRun.write(reinterpret_cast<char *>(runoff), sizeof(float)*numLand);
  fileRun.close();
  delete [] runoff;
  
  //  Discharge from landscape elements
  float * discharge = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"dew");
  strcat(fileName,timeName);
  fileDew.open(fileName,ios::out | ios::binary);
  if (fileDew == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (Dew[k].GetDischarge() > missingData) 
          discharge[k] = (float) (Dew[k].GetDischarge());
      else
          discharge[k] = (float) noData;
      k++;
  }
/*        else {
          discharge[ELEMENT(i,j)] = noData;
          }
          }
          else {
          discharge[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileDew.write(reinterpret_cast<char *>(discharge), sizeof(float)*numLand);
  fileDew.close();
  delete [] discharge;
  
  //  Soil moisture deficit in landscape elements
  float * hsd = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"hsd");
  strcat(fileName,timeName);
  fileHsd.open(fileName,ios::out | ios::binary);
  if (fileHsd == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
    if (Dew[k].GetHbvSoilMoistureDeficit() > missingData) {
      if (Dew[k].GetHbvSoilMoistureDeficit() >= 0.0)
        hsd[k] = (float) ((Dew[k].GetHbvSoilMoistureDeficit())*1000.0);
      else
        hsd[k] = 0.0;
    }
    else
      hsd[k] = (float) noData;
    k++;
  }
/*        else {
          hsd[ELEMENT(i,j)] = noData;
          }
          }
          else {
          hsd[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileHsd.write(reinterpret_cast<char *>(hsd), sizeof(float)*numLand);
  fileHsd.close();
  delete [] hsd;
  
  //  Soil moisture in landscape elements
  float * hsm = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"hsm");
  strcat(fileName,timeName);
  fileHsm.open(fileName,ios::out | ios::binary);
  if (fileHsm == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
    if (Dew[k].GetHbvSoilMoisture() > missingData) {
      if (Dew[k].GetHbvSoilMoisture() >= 0.0)
        hsm[k] = (float) ((Dew[k].GetHbvSoilMoisture())*1000.0);
      else
        hsm[k] = 0.0;
    }
    else
      hsm[k] = (float) noData;
    k++;
  }
/*        else {
          hsm[ELEMENT(i,j)] = noData;
          }
          }
          else {
          hsm[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileHsm.write(reinterpret_cast<char *>(hsm), sizeof(float)*numLand);
  fileHsm.close();
  delete [] hsm;
  
  //  Groundwater in landscape elements
  float * hgw = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"hgw");
  strcat(fileName,timeName);
  fileHgw.open(fileName,ios::out | ios::binary);
  if (fileHgw == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
    //        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (Dew[k].GetHbvUpperZone() > missingData && Dew[k].GetHbvLowerZone() > missingData) 
          hgw[k] = (float) ((Dew[k].GetHbvUpperZone() + Dew[k].GetHbvLowerZone())*1000.0);
      else
          hgw[k] = (float) noData;
      k++;
  }
/*        else {
          hgw[ELEMENT(i,j)] = noData;
          }
          }
          else {
          hgw[ELEMENT(i,j)] = noData;
          }
          }
          }*/
  fileHgw.write(reinterpret_cast<char *>(hgw), sizeof(float)*numLand);
  fileHgw.close();
  delete [] hgw;
  
  //  Lower zone groundwater in landscape elements
  float * hlz = new float [numLand];
  strcpy(fileName,dirName);
  strcat(fileName,"hlz");
  strcat(fileName,timeName);
  fileHlz.open(fileName,ios::out | ios::binary);
  if (fileHlz == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  /*  for (i=0; i<nRows; i++) {
      for (j=0; j<nCols; j++) {*/
  while (k<numLand) {
//        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
      if (Dew[k].GetHbvLowerZone() > missingData) 
          hlz[k] = (float) ((Dew[k].GetHbvLowerZone())*1000.0);
      else
          hlz[k] = (float) noData;
      k++;
  }
/*        else {
          hlz[ELEMENT(i,j)] = noData;
          }
          }
          else {
          hlz[ELEMENT(i,j)] = noData;
          }
          }
  }*/
  fileHlz.write(reinterpret_cast<char *>(hlz), sizeof(float)*numLand);
  fileHlz.close();
  delete [] hlz;
}


void WriteBinaryGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,  
                     int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
  int i,j,k;
  char fileName[100];
  char timeName[100];
  unsigned short int noData16=65535;
  ofstream fileTemp, filePre, fileEva, fileSwe, fileScf, fileSmw, fileRun, fileDew, fileHsd, fileHsm, fileHuz, fileHlz;

  sprintf(timeName,"_%04d_%02d_%02d.bil",datetime.getYear(),datetime.getMonth(),datetime.getDay());

  //  Temperature in landscape elements 
  unsigned short int * temp = new unsigned short int [nRows*nCols];
  strcpy(fileName,"tem");
  strcat(fileName,timeName);
  fileTemp.open(fileName,ios::out | ios::binary);
  if (fileTemp == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedTemperature() > missingData) 
          if (Dew[k].GetTemperature() > missingData) 
            temp[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetTemperature()*10.0)+2731); //temp10K[ELEMENT(i,j)]-2731)/10.0
          else
            temp[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          temp[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        temp[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileTemp.write(reinterpret_cast<char *>(temp), sizeof(unsigned short int)*nRows*nCols);
  fileTemp.close();
  delete [] temp;
  
  //  Precipitation in landscape elements 
  unsigned short int * pre = new unsigned short int [nRows*nCols];
  strcpy(fileName,"pre");
  strcat(fileName,timeName);
  filePre.open(fileName,ios::out | ios::binary);
  if (filePre == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetPrecipitation() > missingData) 
            pre[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetPrecipitation())*1000.0);
          else
            pre[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          pre[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        pre[ELEMENT(i,j)] = noData16;
      }
    }
  }
  filePre.write(reinterpret_cast<char *>(pre), sizeof(unsigned short int)*nRows*nCols);
  filePre.close();
  delete [] pre;
  
  //  Evapotranspiration from landscape elements
  unsigned short int * eva = new unsigned short int [nRows*nCols];
  strcpy(fileName,"eva");
  strcat(fileName,timeName);
  fileEva.open(fileName,ios::out | ios::binary);
  if (fileEva == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetInterceptionLoss() > missingData &&
              Dew[k].GetTranspSoilEvap() > missingData &&
              Dew[k].GetLakeEvap() > missingData)
            eva[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetInterceptionLoss() + Dew[k].GetTranspSoilEvap() + 
                                                       Dew[k].GetLakeEvap())*1000.0);
          else
            eva[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          eva[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        eva[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileEva.write(reinterpret_cast<char *>(eva), sizeof(unsigned short int)*nRows*nCols);
  fileEva.close();
  delete [] eva;
  
  //  Snow store in landscape elements
  unsigned short int * swe = new unsigned short int [nRows*nCols];
  strcpy(fileName,"swe");
  strcat(fileName,timeName);
  fileSwe.open(fileName,ios::out | ios::binary);
  if (fileSwe == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetSnowStore() > missingData) 
            swe[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetSnowStore())*1000.0);
          else
            swe[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          swe[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        swe[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileSwe.write(reinterpret_cast<char *>(swe), sizeof(unsigned short int)*nRows*nCols);
  fileSwe.close();
  delete [] swe;
  
  //  Snowcover fraction in landscape elements
  unsigned short int * scf = new unsigned short int [nRows*nCols];
  strcpy(fileName,"scf");
  strcat(fileName,timeName);
  fileScf.open(fileName,ios::out | ios::binary);
  if (fileScf == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetSnowCoverFraction() > missingData) 
            scf[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetSnowCoverFraction())*1000.0);
          else
            scf[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          scf[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        scf[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileScf.write(reinterpret_cast<char *>(scf), sizeof(unsigned short int)*nRows*nCols);
  fileScf.close();
  delete [] scf;
   
  //  Snow meltwater in landscape elements
  unsigned short int * smw = new unsigned short int [nRows*nCols];
  strcpy(fileName,"smw");
  strcat(fileName,timeName);
  fileSmw.open(fileName,ios::out | ios::binary);
  if (fileSmw == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetMeltWater() > missingData) 
            smw[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetMeltWater())*1000.0);
          else
            smw[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          smw[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        smw[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileSmw.write(reinterpret_cast<char *>(smw), sizeof(unsigned short int)*nRows*nCols);
  fileSmw.close();
  delete [] smw;
  
  //  Runoff from landscape elements
  unsigned short int * runoff = new unsigned short int [nRows*nCols];
  strcpy(fileName,"run");
  strcat(fileName,timeName);
  fileRun.open(fileName,ios::out | ios::binary);
  if (fileRun == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetRunoff() > missingData) 
            runoff[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetRunoff())*1000.0);
          else
            runoff[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          runoff[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        runoff[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileRun.write(reinterpret_cast<char *>(runoff), sizeof(unsigned short int)*nRows*nCols);
  fileRun.close();
  delete [] runoff;
  
  //  Discharge from landscape elements
  unsigned short int * discharge = new unsigned short int [nRows*nCols];
  strcpy(fileName,"dew");
  strcat(fileName,timeName);
  fileDew.open(fileName,ios::out | ios::binary);
  if (fileDew == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetDischarge() > missingData) 
            discharge[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetDischarge())*1000.0);
          else
            discharge[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          discharge[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        discharge[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileDew.write(reinterpret_cast<char *>(discharge), sizeof(unsigned short int)*nRows*nCols);
  fileDew.close();
  delete [] discharge;
  
  //  Soil moisture deficit in landscape elements
  unsigned short int * hsd = new unsigned short int [nRows*nCols];
  strcpy(fileName,"hsd");
  strcat(fileName,timeName);
  fileHsd.open(fileName,ios::out | ios::binary);
  if (fileHsd == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetHbvSoilMoistureDeficit() > missingData) 
            hsd[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetHbvSoilMoistureDeficit())*1000.0);
          else
            hsd[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          hsd[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        hsd[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileHsd.write(reinterpret_cast<char *>(hsd), sizeof(unsigned short int)*nRows*nCols);
  fileHsd.close();
  delete [] hsd;
  
  //  Soil moisture in landscape elements
  unsigned short int * hsm = new unsigned short int [nRows*nCols];
  strcpy(fileName,"hsm");
  strcat(fileName,timeName);
  fileHsm.open(fileName,ios::out | ios::binary);
  if (fileHsm == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetHbvSoilMoisture() > missingData) 
            hsm[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetHbvSoilMoisture())*1000.0);
          else
            hsm[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          hsm[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        hsm[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileHsm.write(reinterpret_cast<char *>(hsm), sizeof(unsigned short int)*nRows*nCols);
  fileHsm.close();
  delete [] hsm;
  
  //  Upper zone groundwater in landscape elements
  unsigned short int * huz = new unsigned short int [nRows*nCols];
  strcpy(fileName,"huz");
  strcat(fileName,timeName);
  fileHuz.open(fileName,ios::out | ios::binary);
  if (fileHuz == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetHbvUpperZone() > missingData) 
            huz[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetHbvUpperZone())*1000.0);
          else
            huz[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          huz[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        huz[ELEMENT(i,j)] = noData16;
      }
    }
    fout << endl;
  }
  fout << endl;
  fileHuz.write(reinterpret_cast<char *>(huz), sizeof(unsigned short int)*nRows*nCols);
  fileHuz.close();
  delete [] huz;
  
  //  Lower zone groundwater in landscape elements
  unsigned short int * hlz = new unsigned short int [nRows*nCols];
  strcpy(fileName,"hlz");
  strcat(fileName,timeName);
  fileHlz.open(fileName,ios::out | ios::binary);
  if (fileHlz == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          //          if (Dew[k].GetAccumulatedPrecipitation() > missingData) 
          if (Dew[k].GetHbvLowerZone() > missingData) 
            hlz[ELEMENT(i,j)] = (unsigned short int) ((Dew[k].GetHbvLowerZone())*1000.0);
          else
            hlz[ELEMENT(i,j)] = noData16;
          k++;
        }
        else {
          hlz[ELEMENT(i,j)] = noData16;
        }
      }
      else {
        hlz[ELEMENT(i,j)] = noData16;
      }
    }
  }
  fileHlz.write(reinterpret_cast<char *>(hlz), sizeof(unsigned short int)*nRows*nCols);
  fileHlz.close();
  delete [] hlz;
}


void WriteAsciiGrid(DistributedElement * const Dew, DateTime datetime, int numLand, int timeStep, int nRows,  
                    int nCols, int noData, double xllCorner, double yllCorner, double cellSize, ofstream &fout)
{
  int i,j,k;
  char fileName[100];
  char dirName[100];
  char timeName[100];

  sprintf(dirName,"./%04d/",datetime.getYear());
  sprintf(timeName,"_%04d_%02d_%02d.asc",datetime.getYear(),datetime.getMonth(),datetime.getDay());

  //  fout << "Temperature in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"tem");
  strcat(fileName,timeName);
  ofstream fileTemp(fileName);
  if (fileTemp == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileTemp << "ncols         " << nCols << endl;
  fileTemp << "nrows         " << nRows << endl;
  fileTemp << "xllcorner     " << xllCorner << endl;
  fileTemp << "yllcorner     " << yllCorner << endl;
  fileTemp << "cellsize      " << cellSize << endl;
  fileTemp << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetTemperature() > missingData) {
            fileTemp.width(15); fileTemp.precision(5); fileTemp.setf(ios::showpoint); fileTemp.setf(ios::fixed); 
            fileTemp << (Dew[k].GetTemperature()) << endl;
          }
          else {
            fileTemp.width(15); fileTemp << noData << endl;
          }
          k++;
        }
        else {
          fileTemp.width(15); fileTemp << noData << endl;
        }
      }
      else {
        fileTemp.width(15); fileTemp << noData << endl;
      }
    }
  }
  fileTemp.close();

  //  fout << "Precipitation in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"pre");
  strcat(fileName,timeName);
  ofstream filePrec(fileName);
  if (filePrec == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  filePrec << "ncols         " << nCols << endl;
  filePrec << "nrows         " << nRows << endl;
  filePrec << "xllcorner     " << xllCorner << endl;
  filePrec << "yllcorner     " << yllCorner << endl;
  filePrec << "cellsize      " << cellSize << endl;
  filePrec << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetPrecipitation() > missingData) {
            filePrec.width(15); filePrec.precision(5); filePrec.setf(ios::showpoint); filePrec.setf(ios::fixed); 
            filePrec << (Dew[k].GetPrecipitation())*1000.0 << endl;
          }
          else {
            filePrec.width(15); filePrec << noData << endl;
          }
          k++;
        }
        else {
          filePrec.width(15); filePrec << noData << endl;
        }
      }
      else {
        filePrec.width(15); filePrec << noData << endl;
      }
    }
  }
  filePrec.close();

  //  fout << "Evapotranspiration from landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"eva");
  strcat(fileName,timeName);
  ofstream fileEvap(fileName);
  if (fileEvap == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileEvap << "ncols         " << nCols << endl;
  fileEvap << "nrows         " << nRows << endl;
  fileEvap << "xllcorner     " << xllCorner << endl;
  fileEvap << "yllcorner     " << yllCorner << endl;
  fileEvap << "cellsize      " << cellSize << endl;
  fileEvap << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetInterceptionLoss() > missingData &&
              Dew[k].GetTranspSoilEvap() > missingData &&
              Dew[k].GetLakeEvap() > missingData) {
            fileEvap.width(15); fileEvap.precision(5); fileEvap.setf(ios::showpoint); fileEvap.setf(ios::fixed); 
            fileEvap << (Dew[k].GetInterceptionLoss() + Dew[k].GetTranspSoilEvap() + 
                         Dew[k].GetLakeEvap())*1000.0 << endl;
          }
          else {
            fileEvap.width(15); fileEvap << noData << endl;
          }
         k++;
        }
        else {
          fileEvap.width(15); fileEvap << noData << endl;
        }
      }
      else {
        fileEvap.width(15); fileEvap << noData << endl;
      }
    }
  }
  fileEvap.close();
  
  //  fout << "Snow store in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"swe");
  strcat(fileName,timeName);
  ofstream fileSwe(fileName);
  if (fileSwe == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileSwe << "ncols         " << nCols << endl;
  fileSwe << "nrows         " << nRows << endl;
  fileSwe << "xllcorner     " << xllCorner << endl;
  fileSwe << "yllcorner     " << yllCorner << endl;
  fileSwe << "cellsize      " << cellSize << endl;
  fileSwe << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetSnowStore() > missingData) {
            fileSwe.width(15); fileSwe.precision(5); fileSwe.setf(ios::showpoint); fileSwe.setf(ios::fixed); 
            fileSwe << (Dew[k].GetSnowStore())*1000.0 << endl;
          }
          else {
            fileSwe.width(15); fileSwe << noData << endl;
          }
          k++;
        }
        else {
          fileSwe.width(15); fileSwe << noData << endl;
        }
      }
      else {
        fileSwe.width(15); fileSwe << noData << endl;
      }
    }
  }
  fileSwe.close();
  
  //  fout << "Snowcover fraction in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"scf");
  strcat(fileName,timeName);
  ofstream fileScf(fileName);
  if (fileScf == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileScf << "ncols         " << nCols << endl;
  fileScf << "nrows         " << nRows << endl;
  fileScf << "xllcorner     " << xllCorner << endl;
  fileScf << "yllcorner     " << yllCorner << endl;
  fileScf << "cellsize      " << cellSize << endl;
  fileScf << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetSnowCoverFraction() > missingData) {
            fileScf.width(15); fileScf.precision(5); fileScf.setf(ios::showpoint); fileScf.setf(ios::fixed); 
            fileScf << (Dew[k].GetSnowCoverFraction())*100.0 << endl;
          }
          else {
            fileScf.width(15); fileScf << noData << endl;
          }
          k++;
        }
        else {
          fileScf.width(15); fileScf << noData << endl;
        }
      }
      else {
        fileScf.width(15); fileScf << noData << endl;
      }
    }
  }
  fileScf.close();
  
  //  fout << "Snow meltwater in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"smw");
  strcat(fileName,timeName);
  ofstream fileSmw(fileName);
  if (fileSmw == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileSmw << "ncols         " << nCols << endl;
  fileSmw << "nrows         " << nRows << endl;
  fileSmw << "xllcorner     " << xllCorner << endl;
  fileSmw << "yllcorner     " << yllCorner << endl;
  fileSmw << "cellsize      " << cellSize << endl;
  fileSmw << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetMeltWater() > missingData) {
            fileSmw.width(15); fileSmw.precision(5); fileSmw.setf(ios::showpoint); fileSmw.setf(ios::fixed); 
            fileSmw << (Dew[k].GetMeltWater())*1000.0 << endl;
          }
          else {
            fileSmw.width(15); fileSmw << noData << endl;
          }
          k++;
        }
        else {
          fileSmw.width(15); fileSmw << noData << endl;
        }
      }
      else {
        fileSmw.width(15); fileSmw << noData << endl;
      }
    }
  }
  fileSmw.close();
  
  //  fout << "Runoff from landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"run");
  strcat(fileName,timeName);
  ofstream fileRunoff(fileName);
  if (fileRunoff == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileRunoff << "ncols         " << nCols << endl;
  fileRunoff << "nrows         " << nRows << endl;
  fileRunoff << "xllcorner     " << xllCorner << endl;
  fileRunoff << "yllcorner     " << yllCorner << endl;
  fileRunoff << "cellsize      " << cellSize << endl;
  fileRunoff << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetRunoff() > missingData) {
            fileRunoff.width(15); fileRunoff.precision(5); fileRunoff.setf(ios::showpoint); fileRunoff.setf(ios::fixed); 
            fileRunoff << (Dew[k].GetRunoff())*1000.0 << endl;
          }
          else {
            fileRunoff.width(15); fileRunoff << noData << endl;
          }
          k++;
        }
        else {
          fileRunoff.width(15); fileRunoff << noData << endl;
        }
      }
      else {
        fileRunoff.width(15); fileRunoff << noData << endl;
      }
    }
  }
  fileRunoff.close();
  
  //  fout << "Discharge from landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"dew");
  strcat(fileName,timeName);
  ofstream fileDisch(fileName);
  if (fileDisch == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileDisch << "ncols         " << nCols << endl;
  fileDisch << "nrows         " << nRows << endl;
  fileDisch << "xllcorner     " << xllCorner << endl;
  fileDisch << "yllcorner     " << yllCorner << endl;
  fileDisch << "cellsize      " << cellSize << endl;
  fileDisch << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetDischarge() > missingData) {
            fileDisch.width(15); fileDisch.precision(5); fileDisch.setf(ios::showpoint); fileDisch.setf(ios::fixed); 
            fileDisch << Dew[k].GetDischarge() << endl;
          }
          else {
            fileDisch.width(15); fileDisch << noData << endl;
          }
          k++;
        }
        else {
          fileDisch.width(15); fileDisch << noData << endl;
        }
      }
      else {
        fileDisch.width(15); fileDisch << noData << endl;
      }
    }
  }
  fileDisch.close();
  
  //  fout << "Soil moisture deficit in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"hsd");
  strcat(fileName,timeName);
  ofstream fileHsd(fileName);
  if (fileHsd == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileHsd << "ncols         " << nCols << endl;
  fileHsd << "nrows         " << nRows << endl;
  fileHsd << "xllcorner     " << xllCorner << endl;
  fileHsd << "yllcorner     " << yllCorner << endl;
  fileHsd << "cellsize      " << cellSize << endl;
  fileHsd << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetHbvSoilMoisture() > missingData) {
            fileHsd.width(15); fileHsd.precision(5); fileHsd.setf(ios::showpoint); fileHsd.setf(ios::fixed);
            fileHsd << (Dew[k].GetHbvSoilMoistureDeficit())*1000.0 << endl;
          }
          else {
            fileHsd.width(15); fileHsd << noData << endl;
          }
          k++;
        }
        else {
          fileHsd.width(15); fileHsd << noData << endl;
        }
      }
      else {
        fileHsd.width(15); fileHsd << noData << endl;
      }
    }
  }
  fileHsd.close();
  
  //  fout << "Soil moisture in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"hsm");
  strcat(fileName,timeName);
  ofstream fileHsm(fileName);
  if (fileHsm == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileHsm << "ncols         " << nCols << endl;
  fileHsm << "nrows         " << nRows << endl;
  fileHsm << "xllcorner     " << xllCorner << endl;
  fileHsm << "yllcorner     " << yllCorner << endl;
  fileHsm << "cellsize      " << cellSize << endl;
  fileHsm << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetHbvSoilMoisture() > missingData) {
            fileHsm.width(15); fileHsm.precision(5); fileHsm.setf(ios::showpoint); fileHsm.setf(ios::fixed);
            fileHsm << (Dew[k].GetHbvSoilMoisture())*1000.0 << endl;
          }
          else {
            fileHsm.width(15); fileHsm << noData << endl;
          }
          k++;
        }
        else {
          fileHsm.width(15); fileHsm << noData << endl;
        }
      }
      else {
        fileHsm.width(15); fileHsm << noData << endl;
      }
    }
  }
  fileHsm.close();
  
  //  fout << "Upper zone in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"huz");
  strcat(fileName,timeName);
  ofstream fileHuz(fileName);
  if (fileHuz == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileHuz << "ncols         " << nCols << endl;
  fileHuz << "nrows         " << nRows << endl;
  fileHuz << "xllcorner     " << xllCorner << endl;
  fileHuz << "yllcorner     " << yllCorner << endl;
  fileHuz << "cellsize      " << cellSize << endl;
  fileHuz << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetHbvUpperZone() > missingData) {
            fileHuz.width(15); fileHuz.precision(5); fileHuz.setf(ios::showpoint); fileHuz.setf(ios::fixed);
            fileHuz << (Dew[k].GetHbvUpperZone())*1000.0 << endl;
          }
          else {
            fileHuz.width(15); fileHuz << noData << endl;
          }
          k++;
        }
        else {
          fileHuz.width(15); fileHuz << noData << endl;
        }
      }
      else {
        fileHuz.width(15); fileHuz << noData << endl;
      }
    }
    fileHuz << endl;
  }
  fileHuz << endl;
  fileHuz.close();
  
  //  fout << "Lower zone in landscape elements at time step " << timeStep << ":\n";
  strcpy(fileName,dirName);
  strcat(fileName,"hlz");
  strcat(fileName,timeName);
  ofstream fileHlz(fileName);
  if (fileHlz == NULL) {
    cout << endl << "Error opening file " << fileName << endl << endl;
    exit (1);
  }
  fileHlz << "ncols         " << nCols << endl;
  fileHlz << "nrows         " << nRows << endl;
  fileHlz << "xllcorner     " << xllCorner << endl;
  fileHlz << "yllcorner     " << yllCorner << endl;
  fileHlz << "cellsize      " << cellSize << endl;
  fileHlz << "NODATA_value  " << noData << endl;
  k=0;
  for (i=0; i<nRows; i++) {
    for (j=0; j<nCols; j++) {
      if (k<numLand) {
        if (Dew[k].GetGeoIndex()==ELEMENT(i,j)) {
          if (Dew[k].GetHbvLowerZone() > missingData) {
            fileHlz.width(15); fileHlz.precision(5); fileHlz.setf(ios::showpoint); fileHlz.setf(ios::fixed);
            fileHlz << (Dew[k].GetHbvLowerZone())*1000.0 << endl;
          }
          else {
            fileHlz.width(15); fileHlz << noData << endl;
          }
          k++;
        }
        else {
          fileHlz.width(15); fileHlz << noData << endl;
        }
      }
      else {
        fileHlz.width(15); fileHlz << noData << endl;
      }
    }
  }
  fileHlz.close();
}


void SetCommonParameters(ParametersCommon * const ParCommonStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  double precGradLow, precGradHigh, gradChangeAltitude, lapseDry, lapseWet; 
  double precCorrRain, precCorrSnow;
  double dayTempMem, lakeEpotPar, kLake, deltaLevel, nLake, maximumLevel;
  double initialSoilMoisture, initialUpperZone, initialLowerZone;
  double saturatedFractionOne, saturatedFractionTwo;
  double initialLakeTemp, initialLakeLevel, initialSnow;
  double initialTotalReservoir;
  int secondsTimestep, numPrec, numTemp, daySnowZero;

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
  finCommonPar.ignore(100,':'); finCommonPar >> secondsTimestep; 
  finCommonPar.ignore(100,':'); finCommonPar >> numPrec; 
  finCommonPar.ignore(100,':'); finCommonPar >> numTemp; 
  finCommonPar.ignore(100,':'); finCommonPar >> precGradLow;
  finCommonPar.ignore(100,':'); finCommonPar >> precGradHigh;
  finCommonPar.ignore(100,':'); finCommonPar >> gradChangeAltitude;
  finCommonPar.ignore(100,':'); finCommonPar >> precCorrRain;
  finCommonPar.ignore(100,':'); finCommonPar >> precCorrSnow;
  finCommonPar.ignore(100,':'); finCommonPar >> lapseDry;
  finCommonPar.ignore(100,':'); finCommonPar >> lapseWet; 
  finCommonPar.ignore(100,':'); finCommonPar >> dayTempMem;
  finCommonPar.ignore(100,':'); finCommonPar >> lakeEpotPar;
  finCommonPar.ignore(100,':'); finCommonPar >> kLake; 
  finCommonPar.ignore(100,':'); finCommonPar >> deltaLevel; 
  finCommonPar.ignore(100,':'); finCommonPar >> nLake; 
  finCommonPar.ignore(100,':'); finCommonPar >> maximumLevel; 
  finCommonPar.ignore(100,':'); finCommonPar >> initialSoilMoisture; 
  finCommonPar.ignore(100,':'); finCommonPar >> initialUpperZone; 
  finCommonPar.ignore(100,':'); finCommonPar >> initialLowerZone; 
  finCommonPar.ignore(100,':'); finCommonPar >> saturatedFractionOne; 
  finCommonPar.ignore(100,':'); finCommonPar >> saturatedFractionTwo; 
  finCommonPar.ignore(100,':'); finCommonPar >> initialLakeTemp; 
  finCommonPar.ignore(100,':'); finCommonPar >> initialLakeLevel; 
  finCommonPar.ignore(100,':'); finCommonPar >> initialSnow; 
  finCommonPar.ignore(100,':'); finCommonPar >> initialTotalReservoir; 
  finCommonPar.ignore(100,':'); finCommonPar >> daySnowZero; 
  ParCommonStore->SetSECONDS_TIMESTEP(secondsTimestep);
  ParCommonStore->SetNUM_PREC_SERIES(numPrec);
  ParCommonStore->SetNUM_TEMP_SERIES(numTemp);
  ParCommonStore->SetPREC_GRAD_LOW(precGradLow);
  ParCommonStore->SetPREC_GRAD_HIGH(precGradHigh);
  ParCommonStore->SetGRAD_CHANGE_ALT(gradChangeAltitude);
  ParCommonStore->SetPREC_CORR_RAIN(precCorrRain);
  ParCommonStore->SetPREC_CORR_SNOW(precCorrSnow);
  ParCommonStore->SetLAPSE_DRY(lapseDry);
  ParCommonStore->SetLAPSE_WET(lapseWet);
  ParCommonStore->SetDAY_TEMP_MEMORY(dayTempMem);
  ParCommonStore->SetLAKE_EPOT_PAR(lakeEpotPar);
  ParCommonStore->SetKLAKE(kLake);
  ParCommonStore->SetDELTA_LEVEL(deltaLevel);
  ParCommonStore->SetNLAKE(nLake);
  ParCommonStore->SetMAXIMUM_LEVEL(maximumLevel);
  ParCommonStore->SetINITIAL_SOIL_MOISTURE(initialSoilMoisture);
  ParCommonStore->SetINITIAL_UPPER_ZONE(initialUpperZone); 
  ParCommonStore->SetINITIAL_LOWER_ZONE(initialLowerZone);
  ParCommonStore->SetINITIAL_SATURATED_ONE(saturatedFractionOne);
  ParCommonStore->SetINITIAL_SATURATED_TWO(saturatedFractionTwo);
  ParCommonStore->SetINITIAL_LAKE_TEMP(initialLakeTemp);
  ParCommonStore->SetINITIAL_LAKE_LEVEL(initialLakeLevel);
  ParCommonStore->SetINITIAL_SNOW(initialSnow);
  ParCommonStore->SetINITIAL_TOTAL_RESERVOIR(initialTotalReservoir);
  ParCommonStore->SetDAY_SNOW_ZERO(daySnowZero);
  //  ParCommonStore->SetNumStations(finCommonPar, numPrec, numTemp); 
  finCommonPar.close();
  fout << "Common parameters: \n";
  fout << ParCommonStore->GetSECONDS_TIMESTEP() << endl;
  fout << ParCommonStore->GetNUM_PREC_SERIES() << endl;
  fout << ParCommonStore->GetNUM_TEMP_SERIES() << endl;
  fout << ParCommonStore->GetPREC_GRAD_LOW() << endl;
  fout << ParCommonStore->GetPREC_GRAD_HIGH() << endl;
  fout << ParCommonStore->GetGRAD_CHANGE_ALT() << endl;
  fout << ParCommonStore->GetPREC_CORR_RAIN() << endl;
  fout << ParCommonStore->GetPREC_CORR_SNOW() << endl;
  fout << ParCommonStore->GetLAPSE_DRY() << endl;
  fout << ParCommonStore->GetLAPSE_WET() << endl;
  fout << ParCommonStore->GetDAY_TEMP_MEMORY() << endl;
  fout << ParCommonStore->GetLAKE_EPOT_PAR() << endl;
  fout << ParCommonStore->GetKLAKE() << endl;
  fout << ParCommonStore->GetDELTA_LEVEL() << endl;
  fout << ParCommonStore->GetNLAKE() << endl;
  fout << ParCommonStore->GetMAXIMUM_LEVEL() << endl;
  fout << ParCommonStore->GetINITIAL_SOIL_MOISTURE() << endl;
  fout << ParCommonStore->GetINITIAL_UPPER_ZONE() << endl;
  fout << ParCommonStore->GetINITIAL_LOWER_ZONE() << endl;
  fout << ParCommonStore->GetINITIAL_SATURATED_ONE() << endl;
  fout << ParCommonStore->GetINITIAL_SATURATED_TWO() << endl;
  fout << ParCommonStore->GetINITIAL_LAKE_TEMP() << endl;
  fout << ParCommonStore->GetINITIAL_LAKE_LEVEL() << endl;
  fout << ParCommonStore->GetINITIAL_SNOW() << endl;
  fout << ParCommonStore->GetINITIAL_TOTAL_RESERVOIR() << endl;
  fout << ParCommonStore->GetDAY_SNOW_ZERO() << endl;
  /*  fout << " Prec.    "; 
      for (i=0; i<ParCommonStore->GetNUM_PREC_SERIES(); i++) {
      fout << ParCommonStore->GetSTATION_ALTITUDE(i) << "  ";
      fout << ParCommonStore->GetSTATION_WEIGHT(i) << "  ";
      }
      fout << endl<< " Temp.    ";
      for (i=0; i<ParCommonStore->GetNUM_TEMP_SERIES(); i++) {
      fout << ParCommonStore->GetSTATION_ALTITUDE(ParCommonStore->GetNUM_PREC_SERIES()+i) << "  ";
      fout << ParCommonStore->GetSTATION_WEIGHT(ParCommonStore->GetNUM_PREC_SERIES()+i) << "  ";
      }*/
  fout << endl;
}


void SetLandSurfaceParameters(ParametersLandSurface * const ParLandSurfaceStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i,j,k;
  double interMax, epotPar, wetPerCorr;
  double accTemp, meltTemp, snowMeltRate, iceMeltRate, freezeEff, maxRel, albedo, cvSnow;

  /*  cout << " File with landsurface parameters: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finLandSurfacePar(fileName);  // Open for reading
  if (finLandSurfacePar == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finLandSurfacePar.getline(buffer, 256);
  for (i=0; i<numberLandSurfaceClasses; i++) {
    finLandSurfacePar >> buffer >> j >> interMax >> epotPar >> wetPerCorr >> accTemp >> meltTemp ;
    finLandSurfacePar >> snowMeltRate >> iceMeltRate >> freezeEff >> maxRel >> albedo >> cvSnow;
    if (i!=j) {
      cout << endl << " Error in land surface parameter file, parameter no. " 
           << i << endl << endl;
      exit (1);
    } 
    ParLandSurfaceStore[i].SetINTER_MAX(interMax);
    ParLandSurfaceStore[i].SetEPOT_PAR(epotPar);
    ParLandSurfaceStore[i].SetWET_PER_CORR(wetPerCorr);
    ParLandSurfaceStore[i].SetACC_TEMP(accTemp);
    ParLandSurfaceStore[i].SetMELT_TEMP(meltTemp);
    ParLandSurfaceStore[i].SetSNOW_MELT_RATE(snowMeltRate);
    ParLandSurfaceStore[i].SetICE_MELT_RATE(iceMeltRate);
    ParLandSurfaceStore[i].SetFREEZE_EFF(freezeEff);
    ParLandSurfaceStore[i].SetMAX_REL(maxRel);
    ParLandSurfaceStore[i].SetALBEDO(albedo);
    ParLandSurfaceStore[i].SetCV_SNOW(cvSnow);
    // Set snow distribution parameters
    SetSnowDistribution(&ParLandSurfaceStore[i], cvSnow);
  }
  finLandSurfacePar.close();

  fout << "Land surface parameters: \n";
  for (i=0; i<numberLandSurfaceClasses; i++) {
    fout << ParLandSurfaceStore[i].GetINTER_MAX() << "    ";
    fout << ParLandSurfaceStore[i].GetEPOT_PAR() << "    ";
    fout << ParLandSurfaceStore[i].GetWET_PER_CORR() << "    ";
    fout << ParLandSurfaceStore[i].GetACC_TEMP() << "    ";
    fout << ParLandSurfaceStore[i].GetMELT_TEMP() << "    ";
    fout << ParLandSurfaceStore[i].GetSNOW_MELT_RATE() << "    ";
    fout << ParLandSurfaceStore[i].GetICE_MELT_RATE() << "    ";
    fout << ParLandSurfaceStore[i].GetFREEZE_EFF() << "    ";
    fout << ParLandSurfaceStore[i].GetMAX_REL() << "    ";
    fout << ParLandSurfaceStore[i].GetALBEDO() << "    ";
    fout << ParLandSurfaceStore[i].GetCV_SNOW() << endl;
    for (k=0; k<numberSnowClasses; k++) {
      fout << ParLandSurfaceStore[i].GetSNOW_WEIGHT(k) << "  ";
    } 
    fout << endl;
  }
  fout << endl;
}


void SetSnowDistribution(ParametersLandSurface * thisParLandSurface, double cvSnow)
{
  int k;
  double stdDevNorm, meanNorm, sumNorm, sumNorm2;
  double stdNormVar[numberSnowClasses], logNormWeight[numberSnowClasses];
  stdNormVar[0]=-2.326347;
  stdNormVar[1]=-1.644853476;
  stdNormVar[2]=-1.036433474;
  stdNormVar[3]=-0.385320604;
  stdNormVar[4]=0.385320604;
  stdNormVar[5]=1.036433474;
  stdNormVar[6]=1.644853476;
  stdNormVar[7]=2.326347;
  stdNormVar[8]=3.71909027;
  //  stdNormVar[8]=4.265043367;
  meanNorm = 0.5*log(1.0/(1.0+cvSnow*cvSnow));
  stdDevNorm = sqrt(log(1.0+cvSnow*cvSnow));
  sumNorm = 0.0;
  for (k=0; k<numberSnowClasses; k++) {
    logNormWeight[k] = exp(stdNormVar[k]*stdDevNorm + meanNorm);
    sumNorm = sumNorm + logNormWeight[k]*probNorm[k];
  } 
  sumNorm2 = 0.0;
  for (k=0; k<numberSnowClasses; k++) {
    logNormWeight[k] = logNormWeight[k]/sumNorm;
    sumNorm2 = sumNorm2 + logNormWeight[k]*probNorm[k];
    thisParLandSurface->SetSNOW_WEIGHT(k,logNormWeight[k]);
    //    cout << k << "  " << logNormWeight[k] << endl;
  } 
  if (sumNorm2 < 1.0-epsilon || sumNorm2 > 1.0+epsilon) {
    cout << endl << " Sum of snow distribution weights = " << sumNorm2 << endl << endl;
    exit(1);
  }
}


void SetSubSurfaceHbvParameters(ParametersSubSurfaceHbv * const ParSubSurfaceHbvStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i, j;
  double fc, fcdel, beta, infmax; 
  double kuz, alfa, perc, klz, draw;

  /*  cout << " File with HBV subsurface parameters: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finSubSurface(fileName);  // Open for reading
  if (finSubSurface == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finSubSurface.getline(buffer, 256);
  for (i=0; i<numberSoilClasses; i++) {
    finSubSurface >> buffer >> j >> fc >> fcdel >> beta >> infmax;
    finSubSurface >> kuz >> alfa >> perc >> klz >> draw;
    if (i!=j) {
      cout << endl << " Error in HBV subsurface parameter file, parameter no. " 
           << i << endl << endl;
      exit (1);
    } 
    ParSubSurfaceHbvStore[i].SetFC(fc);
    ParSubSurfaceHbvStore[i].SetFCDEL(fcdel);
    ParSubSurfaceHbvStore[i].SetBETA(beta);
    ParSubSurfaceHbvStore[i].SetINFMAX(infmax);
    ParSubSurfaceHbvStore[i].SetKUZ(kuz);
    ParSubSurfaceHbvStore[i].SetALFA(alfa);
    ParSubSurfaceHbvStore[i].SetPERC(perc);
    ParSubSurfaceHbvStore[i].SetKLZ(klz);
    ParSubSurfaceHbvStore[i].SetDRAW(draw);
  }
  finSubSurface.close();
  fout << "HBV subsurface parameters: \n";
  for (i=0; i<numberSoilClasses; i++) {
    fout << ParSubSurfaceHbvStore[i].GetFC() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetFCDEL() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetBETA() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetINFMAX() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetKUZ() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetALFA() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetPERC() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetKLZ() << "    ";
    fout << ParSubSurfaceHbvStore[i].GetDRAW() << "\n";
  }
  fout << endl;
}


void SetKiWaParameters(ParametersKiWa * const ParKiWaStore, ifstream &fileControl, ofstream &fout)
{
  char fileName[80];
  char buffer[256];
  int i, j;
  //  double slopeLength, soilDepth, ovPar1, ovPar2; 
  double soilDepth, ovPar1, ovPar2; 
  double tSat0, effPor, kSat0, a, delta, lambdaKw, rootDepth, wiltPoint, eactPar;

  /*  cout << " File with KiWa subsurface parameters: ";
      cin >> fileName;
      cout << endl;*/
  fileControl.ignore(100,':');
  fileControl >> fileName;
  fileControl.ignore(256,'\n');
  ifstream finKiWa(fileName);  // Open for reading
  if (finKiWa == NULL) {
    cout << endl << " Error opening file " << fileName << endl << endl;
    exit(1);
  }
  finKiWa.getline(buffer, 256);
  for (i=0; i<numberSoilClasses; i++) {
    //    finKiWa >> buffer >> j >> slopeLength >> soilDepth >> ovPar1 >> ovPar2 >> tSat0 >> effPor;
    finKiWa >> buffer >> j >> soilDepth >> ovPar1 >> ovPar2 >> tSat0 >> effPor;
    finKiWa >> kSat0 >> a >> delta >> lambdaKw >> rootDepth >> wiltPoint >> eactPar;
    if (i!=j) {
      cout << endl << " Error in KiWa subsurface parameter file, parameter no. " 
           << i << endl << endl;
      exit (1);
    } 
    //    ParKiWaStore[i].SetSLOPE_LENGTH(slopeLength);
    ParKiWaStore[i].SetSOIL_DEPTH(soilDepth);
    ParKiWaStore[i].SetOV_PAR_1(ovPar1);
    ParKiWaStore[i].SetOV_PAR_2(ovPar2);
    ParKiWaStore[i].SetTSAT_0(tSat0);
    ParKiWaStore[i].SetEFF_POR(effPor);
    ParKiWaStore[i].SetKSAT_0(kSat0);
    ParKiWaStore[i].SetA(a);
    ParKiWaStore[i].SetDELTA(delta);
    ParKiWaStore[i].SetLAMBDA_KW(lambdaKw);
    ParKiWaStore[i].SetROOT_DEPTH(rootDepth);
    ParKiWaStore[i].SetWILT_POINT(wiltPoint);
    ParKiWaStore[i].SetEACT_PAR(eactPar);
  }
  finKiWa.close();
  fout <<"KiWa subsurface parameters: \n";
  for (i=0; i<numberSoilClasses; i++) {
    //    fout << ParKiWaStore[i].GetSLOPE_LENGTH() << "    ";
    fout << ParKiWaStore[i].GetSOIL_DEPTH() << "    ";
    fout << ParKiWaStore[i].GetOV_PAR_1() << "    ";
    fout << ParKiWaStore[i].GetOV_PAR_2() << "    ";
    fout << ParKiWaStore[i].GetTSAT_0() << "    ";
    fout << ParKiWaStore[i].GetEFF_POR() << "    ";
    fout << ParKiWaStore[i].GetKSAT_0() << "    ";
    fout << ParKiWaStore[i].GetA() << "    ";
    fout << ParKiWaStore[i].GetDELTA() << "    ";
    fout << ParKiWaStore[i].GetLAMBDA_KW() << "    ";
    fout << ParKiWaStore[i].GetROOT_DEPTH() << "    ";
    fout << ParKiWaStore[i].GetWILT_POINT() << "    ";
    fout << ParKiWaStore[i].GetEACT_PAR() << "\n";
  }
  fout << endl;
}


void Vegetation::WaterBalance(int timeStep)
{
  double timeResolution=1.0;
  precipitation=GetInputElement()->GetInput(0);
  temp = GetInputElement()->GetInput(1);
  //  cout << "    timeStep " << timeStep << "    Vegetation precipitation " << precipitation;
  //  cout << "    Temperature " << temp << endl;

  /* Potential evapotranspiration */
  wetPeriod = 0.0;
  dryPeriod = timeResolution;
  throughFall = 0.0;
  interceptionLoss = 0.0;
  potev = potentialEvap (temp, landSurfacePar->GetEPOT_PAR());

  /* Water input from precipitation and/or snowmelt > 0 or interception remaining from previous time step */
  if (precipitation > 0.0 || prevInterception > 0.0) {
    interceptionStore = prevInterception + precipitation;
    if (potev > 0.0) {
      wetPeriod = interceptionStore / potev;
    }
    else {
      wetPeriod = 0.0;
    }
    //    cout << interceptionStore << " " << potev  << " " << timeResolution  << " " << wetPeriod << " " << dryPeriod << endl;
    if (wetPeriod < timeResolution) {
      interceptionLoss = potev * wetPeriod;
      dryPeriod = timeResolution - (wetPeriod*landSurfacePar->GetWET_PER_CORR());
    }
    else {
      interceptionLoss = potev * timeResolution;
      dryPeriod = 0.0;
    }
    interceptionStore = interceptionStore - interceptionLoss;

    /* If interception store > INTER_MAX, surplus water is infiltrated through the soil surface
       Soil moistured deficit and percolation is calculated */
    //    cout << interceptionStore << " " << potev  << " " << timeResolution  << " " << wetPeriod << " " << dryPeriod << endl;
    if (interceptionStore > landSurfacePar->GetINTER_MAX()) {
      throughFall = interceptionStore - landSurfacePar->GetINTER_MAX();
      interceptionStore = landSurfacePar->GetINTER_MAX();
    }

    /* If interception store < 0 following evaporation at potential rate then dry period > 0, 
       transpiration and soil evaporation will be calculated for dryPeriod */
    if (interceptionStore < 0.0-epsilon*epsilon) {
      printf("    timeStep = %d      interceptionStore = %f\n",timeStep,interceptionStore);
      /*        exit(1);*/
    }
    if (wetPeriod < 0.0) {
      printf("    timeStep = %d      wetPeriod = %f\n",timeStep,wetPeriod);
      /*        exit(1);*/
    }
    if (dryPeriod < 0.0 || dryPeriod > timeResolution) {
      printf("    timeStep = %d      dryPeriod = %f\n",timeStep,dryPeriod);
      /*        exit(1);*/
    }

    /* If 0 <= interception store <= INTER_MAX, no action is necessary
       Infiltration through soil surface = 0, soil moisture deficit and depth of saturated zone is unchanged */
  }
  prevInterception = interceptionStore;
  //  cout << "throughFall " << throughFall << endl;
}


void Snow::WaterBalance(int timeStep, double waterInput)
{
  int i;
  double waterAdded;                                /*  Water added to snow store (m)  */
  double snowMelt=0.0;                              /*  Snowmelt produced in current time step (m)  */
  double reFreeze=0.0;                              /*  Refreeze of meltwater in snow (m)  */
  double initialSwe;                                /*  Snow water equivalent at start of time step (m)  */

  temp = GetInputElement()->GetInput(1);
  //  cout << "timeStep " << timeStep << "             Snow waterInput " << waterInput << "    Temperature " << temp << endl;
  initialSwe = snowStore + meltWater;

  // Uniform snow distribution
  if (landSurfacePar->GetCV_SNOW() == 0.0) {

    /*  Snow accumulation  */
    if (temp < landSurfacePar->GetACC_TEMP() && waterInput > 0.0) { 
      snowStore = snowStore + waterInput;
      waterAdded = 0.0;
    }
    else {
      waterAdded = waterInput;
    }
    
    if (snowStore > 0.0) {
      /*  Water added to meltwater from water input  */
      meltWater = meltWater + waterAdded;

      /*  Refreeze of meltwater in snow  */
      if (temp < landSurfacePar->GetMELT_TEMP() && meltWater > 0.0) {
        reFreeze = landSurfacePar->GetFREEZE_EFF() * landSurfacePar->GetSNOW_MELT_RATE() * (landSurfacePar->GetMELT_TEMP() - temp);
        if (reFreeze > meltWater) {
          reFreeze = meltWater;
        }
        meltWater = meltWater - reFreeze;
        snowStore = snowStore + reFreeze;
      }
      
      /*  Snow store and snow melt  */
      if (temp > landSurfacePar->GetMELT_TEMP() && snowStore > 0.0) {
        snowMelt = landSurfacePar->GetSNOW_MELT_RATE() * (temp - landSurfacePar->GetMELT_TEMP());
        if (snowMelt > snowStore) {
          snowMelt = snowStore;
        }
        snowStore = snowStore - snowMelt;
        meltWater = meltWater + snowMelt;
      }
      
      /* Water output from snow store */
      if (meltWater > landSurfacePar->GetMAX_REL() * snowStore) {
        waterOutput = meltWater - landSurfacePar->GetMAX_REL() * snowStore;
        meltWater = meltWater - waterOutput;
      }
      else {
        waterOutput = 0.0;
      }
    } 
    else {
      waterOutput = waterInput;
    }

    /*  Fraction of area covered by snow  */
    if (snowStore > 0.0)
      snowCoverFraction = 1.0;
    else
      snowCoverFraction = 0.0;     
  } 

  // Lognormal snow distribution
  else  {

    for (i=0; i<numberSnowClasses; i++) {

      /*  Snow accumulation  */
      if (temp < landSurfacePar->GetACC_TEMP() && waterInput > 0.0) { 
        distSnowStore[i] = distSnowStore[i] + waterInput*landSurfacePar->GetSNOW_WEIGHT(i);
        waterAdded = 0.0;
      }
      else {
        waterAdded = waterInput;
      }

      if (distSnowStore[i] > 0.0) {
        /*  Water added to meltwater from water input  */
        distMeltWater[i] = distMeltWater[i] + waterAdded;

        /*  Refreeze of meltwater in snow  */
        if (temp < landSurfacePar->GetMELT_TEMP() && distMeltWater[i] > 0.0) {
          reFreeze = landSurfacePar->GetFREEZE_EFF() * landSurfacePar->GetSNOW_MELT_RATE() * (landSurfacePar->GetMELT_TEMP() - temp);
          if (reFreeze > distMeltWater[i]) {
            reFreeze = distMeltWater[i];
          }
          distMeltWater[i] = distMeltWater[i] - reFreeze;
          distSnowStore[i] = distSnowStore[i] + reFreeze;
        }
        
        /*  Snow store and snow melt  */
        if (temp > landSurfacePar->GetMELT_TEMP() && distSnowStore[i] > 0.0) {
          snowMelt = landSurfacePar->GetSNOW_MELT_RATE() * (temp - landSurfacePar->GetMELT_TEMP());
          if (snowMelt > distSnowStore[i]) {
            snowMelt = distSnowStore[i];
          }
          distSnowStore[i] = distSnowStore[i] - snowMelt;
          distMeltWater[i] = distMeltWater[i] + snowMelt;
        }
        
        /* Water output from snow store */
        if (distMeltWater[i] > landSurfacePar->GetMAX_REL() * distSnowStore[i]) {
          distWaterOutput[i] = distMeltWater[i] - landSurfacePar->GetMAX_REL() * distSnowStore[i];
          distMeltWater[i] = distMeltWater[i] - distWaterOutput[i];
        }
        else {
          distWaterOutput[i] = 0.0;
        }
      } 
      else {
        distWaterOutput[i] = waterInput;
      }
    }
      
    /*  Accumulate values from distributed snow classes  */
    snowStore = 0.0;
    waterOutput = 0.0;
    meltWater = 0.0;
    for (i=0; i<numberSnowClasses; i++) {
      //      cout << " snow.h " << i << "  distSnowStore " << distSnowStore[i] << "  distMeltWater " << distMeltWater[i] << endl;
      snowStore = snowStore + distSnowStore[i]*probNorm[i];
      waterOutput = waterOutput + distWaterOutput[i]*probNorm[i];
      meltWater = meltWater + distMeltWater[i]*probNorm[i];
    }

    /*  Fraction of area covered by snow  */
    snowCoverFraction = 0.0; 
    for (i=numberSnowClasses-1; i>=0; i--) {
      if (distSnowStore[i] > 0.0) snowCoverFraction = snowCoverFraction + probNorm[i];
    }
    
  }

  snowWaterEquivalentChange = snowStore + meltWater - initialSwe;

}


// ** Algorithm to be performed in case: no input to landscape element from upstream elements
//void HBV::WaterBalance(int timeStep, double waterInput, double snowCoverFraction, double dryPeriod)
// ** Algorithm to be performed in case: input to landscape element from upstream elements
void HBV::WaterBalance(int timeStep, double waterInput, double snowCoverFraction, double dryPeriod, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge)
{
  int i, hoursPerTimeStep;
  double inSoilHour, inSoil, outSoil;
  double lowerPercolation, upperPercolation, upperTemporary;
  double upperZoneMax, lowerZoneMax, drawUp, directRunoff;
  temp = GetInputElement()->GetInput(1);
  //  cout << "timeStep " << timeStep << " " << soilMoisture << " " << upperZone << " " << lowerZone << " " << endl;
  //  cout << "\ntimeStep " << timeStep << "   subSurfacePar->GetFC()  " << subSurfacePar->GetFC() << endl;

  // ** Algorithm to be performed in case: input to landscape element from upstream elements through soil surface
  /* Water input from upland elements added to local water input from vegetation and snow store, discharge (m3/s) -> runoff (m) */
  //  cout << timeStep << "  " << waterInput*1000 << "  " << temp << "\n  " ;
  waterInput = waterInput + (upLandAccumulatedLowerDischarge+upLandAccumulatedUpperDischarge)*commonPar->GetSECONDS_TIMESTEP()/GetLandScapeElement()->GetArea();
  // ** End algorithm to be performed in case: input to landscape element from upstream elements through soil surface

  // ** Algorithm to be performed in case: input to landscape element from upstream elements through subsurface
  /* Water input from upland elements added to local water input from vegetation and snow store, discharge (m3/s) -> runoff (m) */
  //  upperZone = upperZone + upLandAccumulatedUpperDischarge*commonPar->GetSECONDS_TIMESTEP()/GetLandScapeElement()->GetArea();
  //  lowerZone = lowerZone + upLandAccumulatedLowerDischarge*commonPar->GetSECONDS_TIMESTEP()/GetLandScapeElement()->GetArea();
  // ** End algorithm to be performed in case: input to landscape element from upstream elements through subsurface

  /* Maximum upper zone storage */
  upperZoneMax = commonPar->GetMAXIMUM_LEVEL();
  //  upperZoneMax = pow(subSurfacePar->GetINFMAX()/subSurfacePar->GetKUZ(),(1/subSurfacePar->GetALFA()));

  /* Maximum lower zone storage */
  lowerZoneMax = subSurfacePar->GetPERC()/subSurfacePar->GetKLZ();

  /* Capillary rise from lower zone to soil moisture zone */
  drawUp = (2.0*subSurfacePar->GetDRAW())*(lowerZone/lowerZoneMax)*(subSurfacePar->GetFC()-soilMoisture)/(subSurfacePar->GetFC());    
  lowerZone = lowerZone - drawUp;
  if (lowerZone < 0.0) {
    drawUp = drawUp + lowerZone;
    lowerZone = 0.0;
  }
  soilMoisture = soilMoisture + drawUp;
    
  /* Water exceeding lower zone maximum to upper zone */
  if (lowerZone > lowerZoneMax) {
    upperZone = upperZone + lowerZone - lowerZoneMax;
    lowerZone = lowerZoneMax;
  }

  /* Water exceeding upper zone maximum to runoff */
  if (upperZone > upperZoneMax) {
    directRunoff = upperZone - upperZoneMax;
    upperZone = upperZoneMax;
  }
  else {
    directRunoff = 0.0;
  }

  /* Soil moisture deficit in use for glacier-free areas */
  if (dryPeriod > missingData) {
    /* Infiltration to soil moisture zone */
    if (waterInput > subSurfacePar->GetINFMAX()) {
      upperPercolation = waterInput - subSurfacePar->GetINFMAX();
      inSoil = subSurfacePar->GetINFMAX();
    }
    else {
      upperPercolation = 0.0;
      inSoil = waterInput;
    }

//    cout << " subSurfaceHbv   snowCoverFraction   " << snowCoverFraction << endl;
    /* Transpiration and soil evaporation for snow free areas */
    if (dryPeriod > 0.0) {
//    if (dryPeriod > 0.0 && GetLandScapeElement()->GetSnowStore() <= 0.0) {
//    if (dryPeriod > 0.0 && GetLandScapeElement()->GetSnowStore() > 0.0) {
      transpSoilEvap = dryPeriod * (1.0 - snowCoverFraction) * 
                       HBVTranspSoilEvap (soilMoisture, temp, landSurfacePar->GetEPOT_PAR(), 
                                          subSurfacePar->GetFC(), subSurfacePar->GetFCDEL());
    }
    else {
      transpSoilEvap = 0.0;
    }
    soilMoisture = soilMoisture - transpSoilEvap;
    if (soilMoisture < 0.0) {
      transpSoilEvap = transpSoilEvap + soilMoisture;
      soilMoisture = 0.0;
    }
    
    /* Water balance for soil moisture zone per hour, percolation to upper groundwater zone */
    if (inSoil > 0.0) {
      hoursPerTimeStep = commonPar->GetSECONDS_TIMESTEP() / minimumTimeStep;
      inSoilHour = inSoil / hoursPerTimeStep;
      //      cout << " hoursPerTimeStep : " << hoursPerTimeStep << " * " << minimumTimeStep << " = " << commonPar->GetSECONDS_TIMESTEP() << endl;
      //      cout << " inSoilHour : " << inSoilHour << " * " << hoursPerTimeStep << " = " << inSoil << endl;
      for (i=0; i<hoursPerTimeStep; i++) {
        if (soilMoisture < subSurfacePar->GetFC())
          outSoil = inSoilHour*pow(soilMoisture/subSurfacePar->GetFC(),subSurfacePar->GetBETA());
        else
          outSoil = inSoilHour;
        soilMoisture = soilMoisture + inSoilHour - outSoil;
        if (soilMoisture > subSurfacePar->GetFC()) {
          outSoil = outSoil + soilMoisture - subSurfacePar->GetFC();
          soilMoisture = subSurfacePar->GetFC();
        }
        upperPercolation = upperPercolation + outSoil;
      }
    }
    else {
      outSoil = 0.0;
    }
    //    cout << "upperPercolation not Glacier  " << upperPercolation << endl;
  }
  
  /* Soil moisture deficit not in use */
  else {
    /*    transpSoilEvap = 0.0;
          upperPercolation = waterInput;
          soilMoisture = subSurfacePar->GetFC();*/
    transpSoilEvap = 0.0;
    upperPercolation = 0.0;
    inSoil = waterInput;

    /* Water balance for soil moisture zone, percolation to upper groundwater zone */
    if (inSoil > 0.0) {
      if (soilMoisture < subSurfacePar->GetFC())
        outSoil = inSoil*power(soilMoisture/subSurfacePar->GetFC(),subSurfacePar->GetBETA());
      else
        outSoil = inSoil;
      soilMoisture = soilMoisture + inSoil - outSoil;
      if (soilMoisture > subSurfacePar->GetFC()) {
        outSoil = outSoil + soilMoisture - subSurfacePar->GetFC();
        soilMoisture = subSurfacePar->GetFC();
      }
      upperPercolation = upperPercolation + outSoil;
    }
    else {
      outSoil = 0.0;
    }
    //    cout << "Glacier:   upperPercolation " << upperPercolation << "  ";
  }

  /* Water balance for upper groundwater zone, percolation to lower groundwater zone */
  lowerPercolation = subSurfacePar->GetPERC();
  //  if (dryPeriod == missingData) cout << " subSurfacePar->GetPERC() " << subSurfacePar->GetPERC() << "  ";
  //  if (dryPeriod == missingData) cout << " lowerPercolation " << lowerPercolation << "  ";
  upperTemporary = upperZone + 0.5*(upperPercolation-lowerPercolation);
  //  if (dryPeriod == missingData) cout << " upperTemporary " << upperTemporary << "  ";
  if (lowerPercolation > upperZone+upperPercolation) {
    lowerPercolation = upperZone + upperPercolation;
    upperRunoff = 0.0;
  }
  else {
    upperRunoff = subSurfacePar->GetKUZ() * pow(upperTemporary,subSurfacePar->GetALFA());
  }
  upperZone = upperZone + upperPercolation - lowerPercolation - upperRunoff;
  if (upperZone < 0.0) {
    upperRunoff = upperRunoff + upperZone;
    if (upperRunoff < 0.0) {
      lowerPercolation = lowerPercolation + upperRunoff;
      upperRunoff = 0.0;
    }
    if (lowerPercolation < 0.0) lowerPercolation = 0.0;
    upperZone = 0.0;
  }
  //  if (dryPeriod == missingData) cout << " lowerPercolation " << lowerPercolation << "  " << endl;

  /* Water balance for lower groundwater zone */
  lowerZone = lowerZone + lowerPercolation;
  if (lowerZone < 0.0) lowerZone = 0.0;
  lowerRunoff = lowerZone * subSurfacePar->GetKLZ();
  lowerZone = lowerZone - lowerRunoff;

  /* Add direct runoff to upper runoff */
  upperRunoff = upperRunoff + directRunoff;

  /* Runoff from both zones */
  runoff = lowerRunoff + upperRunoff;

  /* Percolation from soil moisture zone to upper zone */
  percSoilUpper = upperPercolation;
}


// ** Algorithm to be performed in case: no input to landscape element from upstream elements
//void LakeWaterBalance::WaterBalance(int timeStep)
// ** Algorithm to be performed in case: input to landscape element from upstream elements
void LakeWaterBalance::WaterBalance(int timeStep, double upLandAccumulatedDischarge)
{
  double tempMemory;
  double stage, newStage;
  double kLake;
  double deltaLevel;
  double nLake;
  double maximumLevel;
  double waterInput;
  precipitation=GetInputElement()->GetInput(0);
  temp = GetInputElement()->GetInput(1);
  /*  cout << "    timeStep " << timeStep << "          Lake precipitation " << precipitation;
      cout << "    Temperature " << temp << endl;*/

  // ** Algorithm to be performed in case: no input to landscape element from upstream elements
  //  waterInput = precipitation;
  // ** Algorithm to be performed in case: input to landscape element from upstream elements
  /* Water input from upland elements added to local water input from vegetation and snow store, discharge (m3/s) -> runoff (m) */
  waterInput = precipitation + upLandAccumulatedDischarge*commonPar->GetSECONDS_TIMESTEP()/GetLandScapeElement()->GetArea();

  /* Maximum lake storage */
  maximumLevel = commonPar->GetMAXIMUM_LEVEL();

  /*  Lake temperature and evaporation  */
  tempMemory = commonPar->GetDAY_TEMP_MEMORY() * 86400.0/commonPar->GetSECONDS_TIMESTEP();
  lakeTemp = lakeTemp*(1.0-(1.0/tempMemory))+temp/tempMemory;
  if (lakeTemp > 0.0)
    lakeEvaporation = commonPar->GetLAKE_EPOT_PAR() * lakeTemp;
  else
    lakeEvaporation = 0.0;

  /*  Lake rating curve parameters  */
  if (GetLandScapeElement()->GetLakeNumber() < 0) {
    deltaLevel = commonPar->GetDELTA_LEVEL();
    nLake = commonPar->GetNLAKE();
    kLake = commonPar->GetKLAKE();
  }
  else {
    deltaLevel = 0.0;
    nLake = 1.0;
    kLake = 1.0;
  }

  /*  Initial lake water level  */
  stage = waterLevel + waterInput - lakeEvaporation;

  /* Water exceeding lake maximum level to runoff */
  if (stage > maximumLevel) {
    runoff = stage - maximumLevel;
    stage = maximumLevel;
  }
  else {
    runoff = 0.0;
  }
    
  /*  Lake outflow  */
  if (stage < (-1)*deltaLevel) {
    runoff = 0.0;
  }
  else {
    runoff = runoff + kLake * pow((stage+deltaLevel),nLake);
  }

  /*  Final runoff  */
  newStage = stage - runoff;
  if (runoff > 0.0) {                                  // New test added in order to allow lake evaporation to draw water below -deltaLevel
    if (newStage + deltaLevel < 0.0) {
      runoff = runoff + newStage + deltaLevel;
      if (runoff < 0.0) runoff = 0.0;
      newStage = (-1)*deltaLevel;
    }
  }

  /*  Final lake water level and lake water level change  */
  waterLevelChange = newStage - waterLevel;
  waterLevel = newStage;
  discharge = runoff*(GetLandScapeElement()->GetArea()*GetLandScapeElement()->GetLake()->GetAreaFraction()/100.0)/commonPar->GetSECONDS_TIMESTEP();

}


void GlacierSurface::WaterBalance(int timeStep, int initialTimeSteps)
{
  //  cout << "start GlacierSurfaceWaterBalance " << endl;
  precipitation=GetInputElement()->GetInput(0);
  temp = GetInputElement()->GetInput(1);
  /*  cout << "    timeStep " << timeStep << "       Glacier precipitation " << precipitation;
      cout << "    Temperature " << temp << endl;*/

  /*  Glacier ice melt  */
//  if (temp > landSurfacePar->GetMELT_TEMP() && GetLandScapeElement()->GetGlacier()->GetSnowStore() == 0.0)
  if (temp > landSurfacePar->GetMELT_TEMP())
    iceMelt = landSurfacePar->GetICE_MELT_RATE() * landSurfacePar->GetSNOW_MELT_RATE() * (temp - landSurfacePar->GetMELT_TEMP());
  else
    iceMelt = 0.0;
  //  cout << temp << "  " << landSurfacePar->GetMELT_TEMP() << "  " << landSurfacePar->GetICE_MELT_RATE() << "  " << GetLandScapeElement()->GetSnowStore() << endl;
  //  cout << "  end GlacierSurfaceWaterBalance " << endl;
}


void KinematicWave::WaterBalance(int timeStep, double waterInput, double snowCoverFraction, double dry_period, double upLandAccumulatedLowerDischarge, double upLandAccumulatedUpperDischarge)
{
  int i,j;
  double SLOPE_LENGTH, SLOPE_ANGLE;
  temp = GetInputElement()->GetInput(1);
  SLOPE_LENGTH = GetLandScapeElement()->GetSlopeLength();
  SLOPE_ANGLE = GetLandScapeElement()->GetSlopeAngle()*acos(-1.0)/180.0;  

  // ** Algorithm to be performed in case: input to landscape element from upstream elements
  /* Water input from upland elements added to local water input from vegetation and snow store, discharge (m3/s) -> runoff (m) */
  waterInput = waterInput + (upLandAccumulatedLowerDischarge+upLandAccumulatedUpperDischarge)*commonPar->GetSECONDS_TIMESTEP()/GetLandScapeElement()->GetArea();
  // ** End algorithm to be performed in case: input to landscape element from upstream elements

  if (timeStep == 0) {
    for (i = 1; i <= numberCharacteristic; i++) {
      smdef[i] = 0;
      perc[i] = 0;
      len_coord[i] = ((double)(i-1)/(double)(numberCharacteristic-1)) * SLOPE_LENGTH;
      fixed_length[i] = len_coord[i];
      sat_depth[i] = kiWaPar->GetSOIL_DEPTH() * (commonPar->GetINITIAL_SATURATED_ONE() + 
                    (commonPar->GetINITIAL_SATURATED_TWO() - commonPar->GetINITIAL_SATURATED_ONE()) * 
                    (i-1) / (numberCharacteristic - 1));
      fixed_sat[i] = 0;
      upp_tim[i] = 0;
      upp_dep[i] = 0;
      upp_len[i] = 0;
    }
  }
  actev = 0;
  actev_loss = 0;
  mean_perc = 0;
  lower_runoff = 0;
  upper_runoff = 0;
  //  cout << "timeStep " << timeStep << "    Kwa waterInput " << waterInput << endl;
  //  cout << "timeStep " << timeStep << "    i " << "50" << "    smdef " << smdef[50] << endl;
  if (waterInput > 0) {
    for (i = 1; i <= numberCharacteristic; i++) {
      smdef[i] = smdef[i] + waterInput;
      if (smdef[i] >= 0) {
        perc[i] = smdef[i];
        smdef[i] = 0;
      }
      else {
        perc[i] = 0;
      }
      mean_perc = mean_perc + perc[i];
    }
    mean_perc = mean_perc / numberCharacteristic;
  }
  //  cout << "timeStep " << timeStep << "    mean_perc " << mean_perc << endl;
  //  cout << "timeStep " << timeStep << "    sat_depth[50] " << sat_depth[50] << endl;

  if (mean_perc > 0)     /* Kinematic wave with lateral inflow */
    kinematic_wave_with_lateral_inflow (mean_perc, len_coord, sat_depth, &lower_runoff, &upper_runoff, upp_tim, upp_len, upp_dep, &u_high, 
                                        SLOPE_LENGTH, SLOPE_ANGLE, kiWaPar->GetSOIL_DEPTH(), 
                                        kiWaPar->GetOV_PAR_1(), kiWaPar->GetOV_PAR_2(), kiWaPar->GetEFF_POR(),  
                                        kiWaPar->GetKSAT_0(), kiWaPar->GetA());
  else                   /* Kinematic wave without lateral inflow */
    kinematic_wave_without_lateral_inflow (len_coord, sat_depth, &lower_runoff, &upper_runoff, upp_tim, upp_len, upp_dep, &u_high, 
                                        SLOPE_LENGTH, SLOPE_ANGLE, kiWaPar->GetSOIL_DEPTH(), 
                                        kiWaPar->GetOV_PAR_1(), kiWaPar->GetOV_PAR_2(), kiWaPar->GetEFF_POR(), 
                                        kiWaPar->GetKSAT_0(), kiWaPar->GetA());
  //  cout << "timeStep " << timeStep << "    lower_runoff " << lower_runoff << "    upper_runoff " << upper_runoff << endl;
  /*  Redistribution of saturated depth profile to fixed length coordinates  */
  i = numberCharacteristic;
  j = numberCharacteristic; 
  while (len_coord[numberCharacteristic] <= fixed_length[j]) {
    fixed_sat[j] = sat_depth[numberCharacteristic];
    j--;
  }
  while (j > 0) {
    while (len_coord[i] > fixed_length[j] && i > 1) i--;
    if (len_coord[i] > fixed_length[j]) {
      fixed_sat[j] = sat_depth[i];
    }
    else {
      fixed_sat[j] = sat_depth[i] + (fixed_length[j] - len_coord[i]) * (sat_depth[i+1] - sat_depth[i]) / (len_coord[i+1] - len_coord[i]);
    }
    j--;
  }
  for (i=1; i<= numberCharacteristic; i++) {
    len_coord[i] = fixed_length[i];
    if (fixed_sat[i] > 0)
      sat_depth[i] = fixed_sat[i];
    else
      sat_depth[i] = 0;
  }
  /* If dry period > 0, calculate volumetric water content at the soil surface, actual evapotranspiration, 
     saturated depth and soil moisture deficit along each characteristic curve in lower layer.
     If upper layer detph > 0, the water demanded by evapotranspiration is extracted from this layer first */
  if (dry_period > 0) {
    /* Upper layer */
    if (u_high > 0) {
      evol_upper = dry_period * (1.0 - snowCoverFraction) * potentialEvap (temp, landSurfacePar->GetEPOT_PAR());
      for (j = 1; j <= u_high; j++) {
        upp_dep[j] = upp_dep[j] - evol_upper;
      }
    }
    /* Lower layer */
    for (i = 1; i <= numberCharacteristic; i++) {
      gw_h = (kiWaPar->GetSOIL_DEPTH() - sat_depth[i]) / cos(SLOPE_ANGLE);
      field_capacity = kiWaPar->GetTSAT_0() * exp(kiWaPar->GetLAMBDA_KW() * gw_h);
      teta = field_capacity + smdef[i]/kiWaPar->GetROOT_DEPTH();
      /*        if (gw_h < 0-EPS_KW) {
                printf("    gw_h = %f < 0    i = %d    index = %d \n",gw_h,i,index);
                exit(1);
                }*/
      /*        if (teta < kiWaPar->GetWILT_POINT()-kiWaPar->GetEPS_KW() || teta > kiWaPar->GetTSAT_0()+kiWaPar->GetEPS_KW()) {
                if (teta < 0-kiWaPar->GetEPS_KW() || teta > kiWaPar->GetTSAT_0()+kiWaPar->GetEPS_KW()) { 
                printf("\n\n    kiWaPar->GetWILT_POINT() = %f    teta = %f    gw_h = %f    fc = %f    smdef = %f    i = %d    index = %d \n\n",
                kiWaPar->GetWILT_POINT(),teta,gw_h,field_capacity,smdef[i],i,index);
                exit(1);
                }*/
      if (teta < 0) teta = 0;
      if (teta > kiWaPar->GetTSAT_0()) teta = kiWaPar->GetTSAT_0();
      evaporated_volume = dry_period * (1.0 - snowCoverFraction) * KiWaTranspSoilEvap (teta, temp, landSurfacePar->GetEPOT_PAR(), 
                                                            kiWaPar->GetEACT_PAR(), kiWaPar->GetTSAT_0(), kiWaPar->GetWILT_POINT());
      /* Find index of characteristic curve in upper layer with length coordinate <= length coordinate along characteristic curve in lower layer */
      if (u_high > 0) {
        if (len_coord[i] >= upp_len[u_high] ) {
          j = u_high;
          if (j > 1) {
            while (len_coord[i] >= upp_len[j-1] && j > 1) j--;
          }
          if (upp_dep[j] < 0) 
            evaporated_volume = (-1)*upp_dep[j];
          else
            evaporated_volume = 0;
        }
      }
      /* Saturated depth and soil moisture deficit in lower layer */
      if (evaporated_volume > 0) {
        if (gw_h > kiWaPar->GetROOT_DEPTH()) {
          def_par = 1.0 - exp(kiWaPar->GetDELTA() * (gw_h-kiWaPar->GetROOT_DEPTH()));
          gw_h = gw_h + ((evaporated_volume * (1.0 - def_par)) / kiWaPar->GetEFF_POR());
          smdef[i] = smdef[i] - evaporated_volume * def_par;
        }
        else {
          volume_root = kiWaPar->GetEFF_POR() * (kiWaPar->GetROOT_DEPTH() - gw_h);
          if (volume_root >= evaporated_volume) {
            gw_h = gw_h + (evaporated_volume / kiWaPar->GetEFF_POR());
          }
          else {
            gw_h = kiWaPar->GetROOT_DEPTH()  + ((evaporated_volume  - volume_root) / kiWaPar->GetEFF_POR());
          }
        }
        sat_depth[i] = kiWaPar->GetSOIL_DEPTH() - (gw_h * cos(SLOPE_ANGLE));
      }
      actev = actev + (1.0 - snowCoverFraction) * KiWaTranspSoilEvap (teta, temp, landSurfacePar->GetEPOT_PAR(), 
                                          kiWaPar->GetEACT_PAR(), kiWaPar->GetTSAT_0(), kiWaPar->GetWILT_POINT());
    }
    actev = actev / (numberCharacteristic);
    actev_loss = actev * dry_period;
  } 
  /* Soil moisture deficit must be less or equal to zero */
  for (i = 1; i <= numberCharacteristic; i++) {
    if (smdef[i] > 0 + epsilon) {
      printf("    smdef[i] = %f > 0    i = %d    index = %d \n",smdef[i],i,timeStep);
      /*        exit(1);*/
    }
  }
  /* New value of u_high in upper layer */
  if (u_high > 0 ) {
    j = u_high;
    while (upp_dep[j] <= 0 && j > 0) j--;
    u_high = j;
  }
  /* Control saturated depth profile, soil moisture deficit */
  for (i=1; i<= numberCharacteristic; i++) {
    if (sat_depth[i] < 0) sat_depth[i] = 0;
    if (smdef[i] > 0) smdef[i] = 0;
  }
  /* Combine runoff from upper and lower layer */
  runoff = lower_runoff + upper_runoff;

  /* Evapotranspiration */
  transpSoilEvap = actev_loss;

  /* Output data */
  if (timeStep >= GetDateTimeInfo()->GetInitialTimeSteps()) {
    for (i=0; i<GetSelectedKiWaHillslopeElements()->GetNumberElements(); i++) {
      if (GetLandScapeElement()->GetLandIndex() == GetSelectedKiWaHillslopeElements()->GetKiWaElement(i)) {
        KiWaGroundWaterTable(GetLandScapeElement()->GetLandIndex(), timeStep);
        KiWaSoilMoisture(GetLandScapeElement()->GetLandIndex(), timeStep);
      }
    }
  }
  //  if (GetLandScapeElement()->GetLandIndex() == selectedLandIndex) cout << "evapotranspiration     " 
  //              << (GetLandScapeElement()->GetSoilEvaporation() + GetLandScapeElement()->GetInterceptionLoss())*1000 << endl;
  /*  if (GetLandScapeElement()->GetLandIndex() == selectedLandIndex) {
    cout << "groundwater slope    ";
    cout.width(10);cout.precision(5);
    cout << -(kiWaPar->GetSOIL_DEPTH() - sat_depth[numberCharacteristic/2]) << endl;
    }*/
//  cout << GetLandScapeElement()->GetLandIndex() << "   " << kiWaPar->GetNumberSelectedKiWaHillslopeElements() << endl;
  /*  if (GetLandScapeElement()->GetLandIndex() == selectedLandIndex) {
    //    if (timeStep == 943) {
    if (timeStep == 2159) {
      for (i = 1; i <= numberCharacteristic; i++) {
        gw_h = kiWaPar->GetSOIL_DEPTH() - sat_depth[i];
        field_capacity = kiWaPar->GetTSAT_0() * exp(kiWaPar->GetLAMBDA_KW() * gw_h);
        teta = field_capacity + smdef[i]/kiWaPar->GetROOT_DEPTH();
        cout << len_coord[i] << "    " << -gw_h << "    " << teta << endl;
      }
    }
    }*/
}


void KinematicWave::kinematic_wave_with_lateral_inflow (double inflow, double *len_coord, double *sat_depth, double *lower_runoff, double *upper_runoff, 
                                         double *upp_tim, double *upp_len, double *upp_dep, int *u_high, double SLOPE_LENGTH, double SLOPE_ANGLE, 
                                         double SOIL_DEPTH, double OV_PAR_1, double OV_PAR_2, double EFF_POR, double KSAT_0, double A)
{
  int i;
  int acc_num_new = 0;                   /*  Number of characteristic curves in lower layer that has reached the end of the hillslope during current time step  */
  int num_end = 0;                       /*  Number of characteristic curves in lower or upper layer layer that reaches the downslope end during current time step  */
  int new_high = 0;                      /*  Number of new characteristic curves in upper soil layer  */
  
  double sum_end = 0;                     /*  Accumulated discharge from characteristic curves in lower or upper layer layer that reaches 
                                             the downslope end of hillslope during current time step  */
  double remain_distance;                 /*  Remaining distance to downslope end of hillslope  */
  double remain_time;                     /*  Remaining time to downslope end of hillslope  */
  double downslope_depth;                 /*  Saturated depth at downslope end of hillslope in lower or upper layer  */
  double increment_dist;                  /*  Distance increment between new characteristic curves in lower layer  */
  double increment_depth;                 /*  Saturated depth increment between new characteristic curves in lower layer  */
  double inter_time[numberCharacteristic+1];          /*  Time from start of time step to intersection of characteristic curve in lower layer
                                             with level SOIL_DEPTH  */
  double inter_length[numberCharacteristic+1];        /*  Distance from initial position of characteristic curve in lower layer 
                                             in current time step to intersection with level SOIL_DEPTH  */
  double timeResolution=1.0;

  for (i = 1; i <= numberCharacteristic; i++) {
    inter_time[i] = -999;
    inter_length[i] = -999;
  }

  /* Calculate time and length from start of time step to intersection of characteristic curve in lower layer with level SOIL_DEPTH  */
  i = numberCharacteristic+1;
  while (i > 1 && (*u_high + new_high < numberCharacteristic)) {
    i--;
    if (*u_high > 0 && len_coord[i] >= upp_len[*u_high]) {
      inter_time[i] = 0;
      inter_length[i] = len_coord[i];
    }
    else if (sat_depth[i] >= SOIL_DEPTH) {
      new_high++;
      inter_time[i] = 0;
      inter_length[i] = len_coord[i];
      upp_tim[*u_high + new_high] = 0;
      upp_len[*u_high + new_high] = len_coord[i];
      upp_dep[*u_high + new_high] = 0;
    }
    else {
      inter_time[i] = (SOIL_DEPTH - sat_depth[i]) * EFF_POR / inflow;
      inter_length[i] = len_coord[i] + ((KSAT_0 * sin(SLOPE_ANGLE)) / (A*inflow)) * exp(A*(SOIL_DEPTH-sat_depth[i]))
                                      * (1.0 - exp(((-1)*A*inflow/EFF_POR) * inter_time[i]));
      if (inter_time[i] < timeResolution && inter_length[i] < SLOPE_LENGTH) {
        new_high++;
        upp_tim[*u_high + new_high] = inter_time[i];
        upp_len[*u_high + new_high] = inter_length[i];
        upp_dep[*u_high + new_high] = 0;
      }
      else {
        inter_time[i] = -999;
        inter_length[i] = -999;
      }
    }
  }

  /* New values are assigned to index u_high */
  *u_high = *u_high + new_high;
  if (*u_high > numberCharacteristic) {
    printf("    u_high = %d    new_high = %d\n",*u_high,new_high);
    /*    exit(1);*/
  }

  /* Calculate time of arrival at downslope end for each characteristic curve in lower layer.
     If time of arrival is within current time step, calculate discharge at downslope end and assign 
     the average of these values to upper layer discharge from hillslope for current time step. */
  for (i = numberCharacteristic; i >= 1; i--) {
    if (inter_time[i] >= 0) {
      remain_distance = SLOPE_LENGTH - inter_length[i];
      remain_time = inter_time[i] + EFF_POR * exp(A*(SOIL_DEPTH-SOIL_DEPTH)) * remain_distance / (KSAT_0 * sin(SLOPE_ANGLE));
      downslope_depth = SOIL_DEPTH;
    }
    else {
      remain_distance = SLOPE_LENGTH - len_coord[i];
      remain_time = EFF_POR * (1.0/(A*inflow)) * (A*(SOIL_DEPTH-sat_depth[i]) - 
                                                  log(exp(A*(SOIL_DEPTH-sat_depth[i])) - A*inflow*remain_distance/(KSAT_0 * sin(SLOPE_ANGLE))));
      downslope_depth = sat_depth[i] + inflow * remain_time / EFF_POR;
    }
    if (remain_time <= timeResolution) {
      num_end = num_end + 1;                                /* New characteristic curve in lower layer */
      sum_end = sum_end + KSAT_0 * sin(SLOPE_ANGLE) * (1.0/A) * exp(A*SOIL_DEPTH) * (1 - exp((-1)*A*downslope_depth));
    }
  }
  if (num_end > 0) {
    *lower_runoff = sum_end / num_end;
  }
  else {
    *lower_runoff = KSAT_0 * sin(SLOPE_ANGLE) * (1.0/A) * exp(A*SOIL_DEPTH) * (1 - exp((-1)*A*sat_depth[numberCharacteristic]));
  }
  *lower_runoff = *lower_runoff / SLOPE_LENGTH;                           /* Convert m2/time_step to m/time_step */
  /*  *lower_runoff = *lower_runoff * 1000 / SLOPE_LENGTH;                           * Convert m2/time_step to mm/time_step */
  /*  *lower_runoff = *lower_runoff * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);        * Convert m2/time_step to m3/second */
  
  /* Update number of new characteristic curves in lower layer to be started */
  acc_num_new = num_end;
  if (acc_num_new > numberCharacteristic) {
    printf("    ERROR    acc_num_new = %d    numberCharacteristic = %d\n",acc_num_new, numberCharacteristic);
    /*    exit(1);*/
  }

  /* Calculate length coordinates along all characteristic curves in lower layer that will not reach end of hillslope during current time step  */
  for (i = numberCharacteristic - acc_num_new; i >= 1; i--) {
    if (inter_time[i] >= 0) {
      len_coord[i] = inter_length[i] + KSAT_0 * sin(SLOPE_ANGLE) * (1.0/EFF_POR) * exp(A*(SOIL_DEPTH-SOIL_DEPTH)) * (timeResolution-inter_time[i]);
      sat_depth[i] = SOIL_DEPTH;
    }
    else {
      len_coord[i] = len_coord[i] + KSAT_0 * sin(SLOPE_ANGLE) * (1.0/(A*inflow)) * exp(A*(SOIL_DEPTH-sat_depth[i])) * 
        (1 - exp((-1)*A*inflow*timeResolution/EFF_POR));
      sat_depth[i] = sat_depth[i] + inflow * timeResolution / EFF_POR;
    }
  }

  if (acc_num_new > 0) {
    /* Update indices for each characteristic curve in lower layer that has not reached the downslope end */
    for (i = numberCharacteristic - acc_num_new; i >= 1; i--) {
      len_coord[i+acc_num_new] = len_coord[i];
      sat_depth[i+acc_num_new] = sat_depth[i];
      inter_time[i+acc_num_new] = inter_time[i];
      inter_length[i+acc_num_new] = inter_length[i];
    }

    /* For each characteristic curve in lower layer that has reached the downslope end, start a new characteristic curve at the upslope end */
    increment_dist = len_coord[1] / (acc_num_new+1);
    increment_depth = sat_depth[1] / (acc_num_new+1);
    for (i = 1; i <= acc_num_new; i++) {
      len_coord[i] = increment_dist * i;
      sat_depth[i] = increment_depth * i;
    }
    acc_num_new = 0;
  }

  /* Calculate discharge from lower layer */
  /*  i = numberCharacteristic;
      *lower_runoff = KSAT_0 * sin(SLOPE_ANGLE) * (1.0/A) * exp(A*SOIL_DEPTH) * (1 - exp((-1)*A*sat_depth[i]));
      *  *lower_runoff = *lower_runoff * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step *
      *lower_runoff = *lower_runoff * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */

  if (*u_high > 0) {
    /* Calculate time of arrival at downslope end for each characteristic curve in the upper layer.
       If time of arrival is within current time step, calculate discharge at downslope end and assign 
       the average of these values to upper layer discharge from hillslope for current time step. */
    num_end = 0;
    sum_end = 0;
    for (i = 1; i <= *u_high; i++) {
      remain_distance = SLOPE_LENGTH - upp_len[i];
      remain_time = (1/inflow) * power((inflow*remain_distance/OV_PAR_1 + power(upp_dep[i],OV_PAR_2)),(1/OV_PAR_2)) - (upp_dep[i]/inflow);
      downslope_depth = upp_dep[i] + inflow * remain_time;
      if (remain_time <= timeResolution) {
        num_end = num_end + 1;                                
        sum_end =sum_end + OV_PAR_1 * power(downslope_depth,OV_PAR_2) / SLOPE_LENGTH;                       /* Convert m2/time_step to m/time_step */
        /*      sum_end =sum_end + OV_PAR_1 * power(downslope_depth,OV_PAR_2) * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step */
        /*        sum_end =sum_end + OV_PAR_1 * power(downslope_depth,OV_PAR_2) * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */
      }
    }
    if (num_end > 0) 
      *upper_runoff = sum_end / num_end;
    else
      *upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) / SLOPE_LENGTH;                       /* Convert m2/time_step to m/time_step */
    /*      *upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step */
      /*      *upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */
    
    /* Update indices for each characteristic curve in the upper layer that has not reached the downslope end and index u_high */
    i = 1;
    while (i+num_end <= *u_high) {
      upp_tim[i] = upp_tim[i + num_end];
      upp_len[i] = upp_len[i + num_end];
      upp_dep[i] = upp_dep[i + num_end];
      i++;
    }
    *u_high = *u_high - num_end;
    if (*u_high < 0) *u_high = 0;
    
    /* Calculate length coordinates along all characteristic curves in the upper layer that 
       will not reach the end of the hillslope during the current time step  */
    for (i = 1; i <= *u_high; i++) {
      upp_len[i] = upp_len[i] + (OV_PAR_1 / (inflow)) * (power((upp_dep[i] + inflow*(timeResolution-upp_tim[i])),OV_PAR_2) - power(upp_dep[i],OV_PAR_2));
      upp_dep[i] = upp_dep[i] + inflow * (timeResolution-upp_tim[i]);
      upp_tim[i] = 0;
      if (upp_len[i] > SLOPE_LENGTH) {
        /*      printf("\n    upp_len[i] = %f    i = %d\n\n",upp_len[i],i);*/
        upp_len[i] = SLOPE_LENGTH;
      }
    }
  }

  return;
}


void KinematicWave::kinematic_wave_without_lateral_inflow (double *len_coord, double *sat_depth, double *lower_runoff, double *upper_runoff, 
                                            double *upp_tim, double *upp_len, double *upp_dep, int *u_high, double SLOPE_LENGTH, double SLOPE_ANGLE, 
                                            double SOIL_DEPTH, double OV_PAR_1, double OV_PAR_2, double EFF_POR, double KSAT_0, double A)
{
  int i;
  int acc_num_new = 0;                   /*  Number of characteristic curves in lower layer that has reached the end of the hillslope during current time step  */
  int num_end = 0;                       /*  Number of characteristic curves in lower or upper layer that reaches the downslope end during current time step  */

  double sum_end = 0;                     /*  Accumulated discharge from characteristic curves in lower or upper layer layer that reaches 
                                              the downslope end of hillslope during current time step  */
  double remain_distance;                 /*  Remaining distance to downslope end of hillslope  */
  double remain_time;                     /*  Remaining time to downslope end of hillslope  */
  double increment_dist;                  /*  Distance increment between new characteristic curves in lower layer  */
  double increment_depth;                 /*  Saturated depth increment between new characteristic curves in lower layer  */
  double timeResolution=1.0;

  /* Calculate time of arrival at downslope end for each characteristic curve in lower layer.
     If time of arrival is within current time step, calculate discharge at downslope end and assign 
     the average of these values to upper layer discharge from hillslope for current time step. */
  for (i = numberCharacteristic; i >= 1; i--) {
    remain_distance = SLOPE_LENGTH - len_coord[i];
    remain_time = EFF_POR * exp(A*(sat_depth[i]-SOIL_DEPTH)) * remain_distance / (KSAT_0 * sin(SLOPE_ANGLE));
    if (remain_time <= timeResolution) {
      num_end = num_end + 1;                                /* New characteristic curve in lower layer */
      sum_end = sum_end + KSAT_0 * sin(SLOPE_ANGLE) * (1.0/A) * exp(A*SOIL_DEPTH) * (1 - exp((-1)*A*sat_depth[i]));
    }
  }
  if (num_end > 0) {
    *lower_runoff = sum_end / num_end;
  }
  else {
    *lower_runoff = KSAT_0 * sin(SLOPE_ANGLE) * (1.0/A) * exp(A*SOIL_DEPTH) * (1 - exp((-1)*A*sat_depth[numberCharacteristic]));
  }
  *lower_runoff = *lower_runoff / SLOPE_LENGTH;                           /* Convert m2/time_step to m/time_step */
  /*  *lower_runoff = *lower_runoff * 1000 / SLOPE_LENGTH;                           * Convert m2/time_step to mm/time_step */
  /*  *lower_runoff = *lower_runoff * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);        * Convert m2/time_step to m3/second */
  
  /* Update number of new characteristic curves in lower layer to be started */
  acc_num_new = num_end;
  if (acc_num_new > numberCharacteristic) {
    printf("    ERROR    acc_num_new = %d    numberCharacteristic = %d\n",acc_num_new, numberCharacteristic);
    /*    exit(1);*/
  }

  /* Calculate length coordinates along all characteristic curves in lower layer that will not reach end of hillslope during current time step.
     Saturated thickness remains unchanged in the absence of lateral inflow */
  for (i = numberCharacteristic - acc_num_new; i >= 1; i--) {
    len_coord[i] = len_coord[i] + KSAT_0 * sin(SLOPE_ANGLE) * (1.0/EFF_POR) * exp(A*(SOIL_DEPTH-sat_depth[i])) * timeResolution;
  }

  if (acc_num_new > 0) {
    /* Update indices for each characteristic curve in lower layer that has not reached the downslope end */
    for (i = numberCharacteristic - acc_num_new; i >= 1; i--) {
      len_coord[i+acc_num_new] = len_coord[i];
      sat_depth[i+acc_num_new] = sat_depth[i];
    }

    /* For each characteristic curve in lower layer that has reached the downslope end, start a new characteristic curve at the upslope end */
    increment_dist = len_coord[1] / (acc_num_new+1);
    increment_depth = sat_depth[1] / (acc_num_new+1);
    for (i = 1; i <= acc_num_new; i++) {
      len_coord[i] = increment_dist * i;
      sat_depth[i] = increment_depth * i;
    }
    acc_num_new = 0;
  }

  /* Calculate discharge from lower layer */
  /*  i = numberCharacteristic;
      *lower_runoff = KSAT_0 * sin(SLOPE_ANGLE) * (1.0/A) * exp(A*SOIL_DEPTH) * (1 - exp((-1)*A*sat_depth[i]));
      *  *lower_runoff = *lower_runoff * 1000 / SLOPE_LENGTH;                           * Convert m2/time_step to mm/time_step *
      *lower_runoff = *lower_runoff * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);        * Convert m2/time_step to m3/second */

  if (*u_high > 0) {
    /* Calculate time of arrival at downslope end for each characteristic curve in the upper layer.
       If time of arrival is within current time step, calculate discharge at downslope end and assign 
       the average of these values to upper layer discharge from hillslope for current time step. */
    num_end = 0;
    sum_end = 0;
    for (i = 1; i <= *u_high; i++) {
      remain_distance = SLOPE_LENGTH - upp_len[i];
      remain_time = remain_distance / (OV_PAR_1 * OV_PAR_2 * power(upp_dep[i],OV_PAR_2-1.0));
      if (remain_time <= timeResolution) {
        num_end = num_end + 1;                                
        sum_end =sum_end + OV_PAR_1 * power(upp_dep[i],OV_PAR_2) / SLOPE_LENGTH;                       /* Convert m2/time_step to m/time_step */
        /*      sum_end =sum_end + OV_PAR_1 * power(upp_dep[i],OV_PAR_2) * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step */
        /*        sum_end =sum_end + OV_PAR_1 * power(upp_dep[i],OV_PAR_2) * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */
      }
    }
    if (num_end > 0) 
      *upper_runoff = sum_end / num_end;
    else
      *upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) / SLOPE_LENGTH;                       /* Convert m2/time_step to m/time_step */
    /*      *upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) * 1000 / SLOPE_LENGTH;                       * Convert m2/time_step to mm/time_step */
      /*      upper_runoff = OV_PAR_1 * power(upp_dep[1],OV_PAR_2) * CATCHMENT_AREA / SLOPE_LENGTH / (3600.0 * timeResolution);    * Convert m2/time_step to m3/second */
    
    /* Update indices for each characteristic curve in the upper layer that has not reached the downslope end and index u_high */
    i = 1;
    while (i+num_end <= *u_high) {
      upp_tim[i] = upp_tim[i + num_end];
      upp_len[i] = upp_len[i + num_end];
      upp_dep[i] = upp_dep[i + num_end];
      i++;
    }
    *u_high = *u_high - num_end;
    if (*u_high < 0) *u_high = 0;
    
    /* Calculate length coordinates along all characteristic curves in the upper layer that 
       will not reach the end of the hillslope during the current time step 
       Depth of upper layer remains unchanged in the absence of lateral inflow */
    for (i = 1; i <= *u_high; i++) {
      upp_tim[i] = 0;
      upp_len[i] = upp_len[i] + OV_PAR_1 * OV_PAR_2 * power(upp_dep[i],OV_PAR_2-1.0) * timeResolution;
      if (upp_len[i] > SLOPE_LENGTH) {
        /*      printf("\n    upp_len[i] = %f    i = %d\n\n",upp_len[i],i);*/
        upp_len[i] = SLOPE_LENGTH;
      }
    }
  }

  return;
}


/* Plot groundwater table depth along the hillslope at a specific time step */
void KinematicWave::KiWaGroundWaterTable(int elementId, int timeStep)
{
  FILE *fp_tmp;
  char fileName[100];
  char timeName[100];
  int i;
  int initialTimeSteps;
  int numberTimeSteps;
  double SLOPE_LENGTH, SLOPE_ANGLE;
  DateTime startModelTime;
  DateTime startSimulationTime;
  DateTime endSimulationTime;
  DateTime datetime;

  SLOPE_LENGTH = GetLandScapeElement()->GetSlopeLength();
  SLOPE_ANGLE = GetLandScapeElement()->GetSlopeAngle()*acos(-1.0)/180.0;  

  initialTimeSteps = GetDateTimeInfo()->GetInitialTimeSteps();
  numberTimeSteps = GetDateTimeInfo()->GetNumberTimeSteps();
  startModelTime = GetDateTimeInfo()->GetStartModelTime();
  startSimulationTime =GetDateTimeInfo()->GetStartSimulationTime();
  endSimulationTime = GetDateTimeInfo()->GetEndSimulationTime();
  datetime = startModelTime + timeStep*GetCommonPar()->GetSECONDS_TIMESTEP();

  sprintf(timeName,"%d_%04d%02d%02d_%02d%02d",elementId,datetime.getYear(),datetime.getMonth(),datetime.getDay(),
          datetime.getHour(),datetime.getMinute());
//  strcpy(fileName,"KiWa_groundwater_");
  strcpy(fileName,"KiWa_hillslope_groundwater_");
  strcat(fileName,timeName);
  strcat(fileName,".txt");

  if ((fp_tmp = fopen(fileName, "w")) == NULL) {
    printf("\n    File %s not found!\n\n",fileName);
    exit(1);
  }
  fprintf(fp_tmp,"# Column 1: KiWa groundwater table depth (m)  %s  time step %d\n",timeName,timeStep);
  for (i = 1; i <= numberCharacteristic; i++) {
    fprintf(fp_tmp," %f    %f\n",len_coord[i],(sat_depth[i]-kiWaPar->GetSOIL_DEPTH()));
  }
  fclose(fp_tmp);

  return;
}


/* Plot volumetric soil moisture content along the hillslope at a specific time step */
void KinematicWave::KiWaSoilMoisture(int elementId, int timeStep)
{
  FILE *fp_tmp;
  char fileName[100];
  char timeName[100];
  int i;
  int initialTimeSteps;
  int numberTimeSteps;
  double gw_h, field_capacity, teta;
  double SLOPE_LENGTH, SLOPE_ANGLE;
  DateTime startModelTime;
  DateTime startSimulationTime;
  DateTime endSimulationTime;
  DateTime datetime;

  SLOPE_LENGTH = GetLandScapeElement()->GetSlopeLength();
  SLOPE_ANGLE = GetLandScapeElement()->GetSlopeAngle()*acos(-1.0)/180.0;  

  initialTimeSteps = GetDateTimeInfo()->GetInitialTimeSteps();
  numberTimeSteps = GetDateTimeInfo()->GetNumberTimeSteps();
  startModelTime = GetDateTimeInfo()->GetStartModelTime();
  startSimulationTime =GetDateTimeInfo()->GetStartSimulationTime();
  endSimulationTime = GetDateTimeInfo()->GetEndSimulationTime();
  datetime = startModelTime + timeStep*GetCommonPar()->GetSECONDS_TIMESTEP();

  sprintf(timeName,"%d_%04d%02d%02d_%02d%02d",elementId,datetime.getYear(),datetime.getMonth(),datetime.getDay(),
          datetime.getHour(),datetime.getMinute());
//  strcpy(fileName,"KiWa_soilmoisture_");
  strcpy(fileName,"KiWa_hillslope_soilmoisture_");
  strcat(fileName,timeName);
  strcat(fileName,".txt");

  if ((fp_tmp = fopen(fileName, "w")) == NULL) {
    printf("\n    File %s not found!\n\n",fileName);
    exit(1);
  }
  fprintf(fp_tmp,"# Column 1: KiWa volumetric soil moisture  %s  time step %d\n",timeName,timeStep);
  for (i = 1; i <= numberCharacteristic; i++) {
    gw_h = (kiWaPar->GetSOIL_DEPTH() - sat_depth[i]);
    field_capacity = kiWaPar->GetTSAT_0() * exp(kiWaPar->GetLAMBDA_KW() * gw_h);
    teta = field_capacity + smdef[i]/kiWaPar->GetROOT_DEPTH();
    fprintf(fp_tmp," %f    %f\n",len_coord[i],teta);
  }
  fclose(fp_tmp);

  return;
}


/* Groundwater table depth at a fraction of total hillslope length (from top to bottom) */
double KinematicWave::GetGroundWaterDepth(double lengthFraction) const
{
  return -(kiWaPar->GetSOIL_DEPTH() - sat_depth[(int)(numberCharacteristic*lengthFraction)]);
}


/* Volumetric soil moisture content at a fraction of total hillslope length (from top to bottom) */
double KinematicWave::GetSoilMoisture(double lengthFraction) const
{
  double gw_h, field_capacity, teta;
  gw_h = (kiWaPar->GetSOIL_DEPTH() - sat_depth[(int)(numberCharacteristic*lengthFraction)]);
  field_capacity = kiWaPar->GetTSAT_0() * exp(kiWaPar->GetLAMBDA_KW() * gw_h);
  teta = field_capacity + smdef[(int)(numberCharacteristic*lengthFraction)]/kiWaPar->GetROOT_DEPTH();
  return teta;
}


// Power function 
double power(double base, double exponent)
{
  if (base > 0.0) 
    return pow(base,exponent);
  else if (base == 0.0) 
    return 0.0;
  else if (base >= -epsilon)
    return 0.0;
  else {
    printf("\n\n    *****     function power     base = %f     exponent = %f     *****\n\n",base,exponent);
    return 0.0;
  }
}


// Leap year
int leapYear(int year)
{
  if (year%400==0) 
    return 1;
  else if (year%100==0) 
    return 0;
  else if (year%4==0) 
    return 1;
  else
    return 0;
}


// Find day number
int dayNumber(int year, int month, int day)
{
  int accDays[] = {0,31,59,90,120,151,181,212,243,273,304,334};
  if (month<3) 
    return accDays[month-1]+day;
  else
    return accDays[month-1]+day+leapYear(year);
}


// Day number to date
void dayNo2Date(int dayNo, int year, int * month, int * day)
{      
  int i;
  int accDays[] = {0,31,59,90,120,151,181,212,243,273,304,334};
  if (dayNo>365+leapYear(year)) dayNo=dayNo-365-leapYear(year);
  if (leapYear(year)) 
    for (i=2; i<12; i++) accDays[i]++;
  i=12;
  while (dayNo<=accDays[i-1]) i--;
  *month=i;
  *day=dayNo-accDays[i-1];
}


// Potential evapotranspiration
double potentialEvap (double temp, double epotPar)
{
  if (temp > 0.0)
    return epotPar * temp;
  else
    return 0.0;
}


// Actual evapotranspiration HBV
double HBVTranspSoilEvap(double soilMoist, double temp, double epotPar, double fieldCapacity, double fcDel)
{
  double epot;
  epot = potentialEvap(temp,epotPar);            
  if (soilMoist > fcDel * fieldCapacity) {
    return epot;
  }
  else if (soilMoist > 0.0) {
    return epot * soilMoist / (fieldCapacity * fcDel);
  }
  else {
    return 0.0;
  }
}


// Actual evapotranspiration KinematicWave
double KiWaTranspSoilEvap(double teta, double temp, double epotPar, double eactPar, double tSat0, double wiltPoint)
{
  /* Condition : eactPar * (tSat0-wiltPoint) - wiltPoint > 0 */
  double epot;
  epot = potentialEvap(temp,epotPar);            
  if (teta > eactPar * (tSat0-wiltPoint)) {
    return epot;
  }
  else if (teta > wiltPoint) {
    return epot * (teta - wiltPoint) / (eactPar * (tSat0-wiltPoint) - wiltPoint);
  }
  else {
    return 0.0;
  }
}
