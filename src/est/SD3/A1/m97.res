Sat Oct 23 22:00:49 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/SD3/A1/dat97.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2

$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr
$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m97.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1251.91122188737        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4568E+02 -1.6372E+01  4.6309E+01 -1.5649E+01  9.9622E+01  1.9799E+01 -2.4115E+01 -1.6525E+02 -4.2732E+01 -3.3339E+01
            -1.3204E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1864.27601594427        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.5690E-01  8.3480E-01  1.0738E+00  1.1070E+00  9.2873E-01  9.5762E-01  1.2271E+00  6.9716E-01  1.1655E+00  1.0676E+00
             1.6753E+00
 PARAMETER:  5.5949E-02 -8.0566E-02  1.7120E-01  2.0162E-01  2.6067E-02  5.6693E-02  3.0465E-01 -2.6073E-01  2.5319E-01  1.6540E-01
             6.1601E-01
 GRADIENT:   2.1123E+02 -1.3074E+01 -1.6697E+01  2.6163E+01  2.8754E+01 -5.7993E+00  6.6827E+00  9.8227E+00  4.3779E+01  7.8519E+00
             9.1168E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1880.13527567312        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      178
 NPARAMETR:  9.3512E-01  5.3397E-01  7.2726E-01  1.2857E+00  6.2601E-01  1.0602E+00  1.9231E+00  1.6409E-01  7.6696E-01  8.7036E-01
             1.5954E+00
 PARAMETER:  3.2925E-02 -5.2741E-01 -2.1847E-01  3.5128E-01 -3.6839E-01  1.5843E-01  7.5396E-01 -1.7073E+00 -1.6532E-01 -3.8847E-02
             5.6713E-01
 GRADIENT:   1.4088E+01  1.3338E+01 -3.4998E+01  1.2829E+01  1.3060E+01  1.3547E+01  3.7333E-01  1.1480E+00 -4.3814E+01  2.5155E+01
             6.1088E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1890.60914243926        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      355
 NPARAMETR:  9.1240E-01  4.1733E-01  7.6935E-01  1.3465E+00  6.0220E-01  1.0120E+00  2.1674E+00  1.5416E-01  8.4975E-01  7.8655E-01
             1.4578E+00
 PARAMETER:  8.3188E-03 -7.7389E-01 -1.6221E-01  3.9750E-01 -4.0717E-01  1.1193E-01  8.7355E-01 -1.7698E+00 -6.2812E-02 -1.4010E-01
             4.7691E-01
 GRADIENT:  -2.8075E+01  1.5939E+01  2.1049E+01 -4.8493E+00 -3.8201E+01 -3.8173E+00  3.1329E-01  9.3084E-01 -1.3769E+01  2.9599E+00
            -1.4519E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1895.28175554971        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  9.1658E-01  1.7658E-01  7.9906E-01  1.4837E+00  5.6285E-01  1.0336E+00  3.2426E+00  8.3543E-02  8.7164E-01  7.7099E-01
             1.5113E+00
 PARAMETER:  1.2896E-02 -1.6340E+00 -1.2432E-01  4.9454E-01 -4.7475E-01  1.3308E-01  1.2764E+00 -2.3824E+00 -3.7379E-02 -1.6008E-01
             5.1294E-01
 GRADIENT:  -6.2824E+00  5.1091E+00  3.6228E+01  1.1779E+01 -5.7200E+01  6.2265E+00 -3.3551E+00  2.0221E-01  5.2443E+00 -1.3267E+00
             1.2227E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1897.09539093124        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.1580E-01  7.4220E-02  7.6763E-01  1.5007E+00  5.4077E-01  1.0279E+00  4.8242E+00  4.1195E-02  8.5071E-01  7.5426E-01
             1.5174E+00
 PARAMETER:  1.2041E-02 -2.5007E+00 -1.6445E-01  5.0590E-01 -5.1476E-01  1.2751E-01  1.6737E+00 -3.0894E+00 -6.1683E-02 -1.8202E-01
             5.1702E-01
 GRADIENT:  -1.3247E+00 -1.1146E+00  2.8381E+00 -3.5035E+01 -9.4435E+00  5.0538E+00 -2.0930E+00  4.8878E-02  5.4617E+00 -7.7662E-01
             1.3330E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1898.31899478906        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  9.1741E-01  4.2657E-02  8.2621E-01  1.5312E+00  5.6780E-01  1.0211E+00  6.0823E+00  4.6498E-02  8.3456E-01  7.7973E-01
             1.4973E+00
 PARAMETER:  1.3795E-02 -3.0546E+00 -9.0904E-02  5.2608E-01 -4.6599E-01  1.2089E-01  1.9054E+00 -2.9683E+00 -8.0847E-02 -1.4881E-01
             5.0363E-01
 GRADIENT:   4.1435E+00 -3.5344E+00  3.3594E+00 -2.2537E+01 -2.9953E+00  2.5312E+00 -6.7388E+00  6.1673E-02  3.9367E+00  8.1613E-01
             5.8683E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1898.60045936221        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1058
 NPARAMETR:  9.1560E-01  4.8108E-02  8.5868E-01  1.5483E+00  5.8477E-01  1.0163E+00  6.0542E+00  7.4614E-02  8.2400E-01  8.0815E-01
             1.4903E+00
 PARAMETER:  1.1819E-02 -2.9343E+00 -5.2359E-02  5.3717E-01 -4.3654E-01  1.1612E-01  1.9007E+00 -2.4954E+00 -9.3581E-02 -1.1301E-01
             4.9900E-01
 GRADIENT:  -2.6335E-01  1.1309E-01 -1.6544E+00 -4.9122E+00  2.0224E+00  4.7879E-01  7.6726E-01  1.4061E-01  2.1787E-01 -6.8002E-01
             7.2098E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1899.51839349747        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1219
 NPARAMETR:  9.1636E-01  4.6525E-02  8.6160E-01  1.5505E+00  5.8338E-01  1.0120E+00  5.8300E+00  1.0000E-02  8.3065E-01  8.1160E-01
             1.4856E+00
 PARAMETER:  1.2651E-02 -2.9678E+00 -4.8967E-02  5.3855E-01 -4.3892E-01  1.1196E-01  1.8630E+00 -7.1306E+00 -8.5552E-02 -1.0875E-01
             4.9580E-01
 GRADIENT:   1.4441E+00 -1.0524E+00  6.4907E+00 -1.1322E+00 -9.5092E+00 -1.3611E+00 -4.0498E+00  0.0000E+00 -5.4696E-01  2.5641E+00
            -1.4078E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1899.59112590146        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1399
 NPARAMETR:  9.1540E-01  4.5649E-02  8.5697E-01  1.5512E+00  5.8284E-01  1.0151E+00  5.8807E+00  1.0000E-02  8.3503E-01  8.0257E-01
             1.4875E+00
 PARAMETER:  1.1605E-02 -2.9868E+00 -5.4357E-02  5.3904E-01 -4.3985E-01  1.1497E-01  1.8717E+00 -6.3115E+00 -8.0284E-02 -1.1994E-01
             4.9712E-01
 GRADIENT:  -4.4062E-01  1.4961E+00 -1.0715E+00 -3.8285E+00  1.9031E+00 -3.5201E-01  1.7454E+00  0.0000E+00 -1.7115E+00 -8.8885E-01
            -2.1893E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1899.71247043340        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1579
 NPARAMETR:  9.1594E-01  4.3021E-02  8.5023E-01  1.5527E+00  5.7904E-01  1.0157E+00  6.0432E+00  2.0735E-02  8.3514E-01  7.9757E-01
             1.4878E+00
 PARAMETER:  1.2193E-02 -3.0461E+00 -6.2249E-02  5.3998E-01 -4.4638E-01  1.1556E-01  1.8989E+00 -3.7760E+00 -8.0155E-02 -1.2619E-01
             4.9730E-01
 GRADIENT:   1.0063E+00  1.7685E+00 -1.9167E+00 -1.0666E+00  1.7921E+00 -1.1764E-01  2.6452E+00  1.0800E-02 -1.8913E+00 -1.1771E+00
            -2.8478E-01

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1899.84852250528        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1719
 NPARAMETR:  9.1581E-01  3.9144E-02  8.5660E-01  1.5530E+00  5.8049E-01  1.0163E+00  6.2610E+00  1.0000E-02  8.2603E-01  7.9331E-01
             1.4918E+00
 PARAMETER:  1.2055E-02 -3.1386E+00 -5.3786E-02  5.3982E-01 -4.4360E-01  1.1630E-01  1.9353E+00 -5.4915E+00 -9.1851E-02 -1.3243E-01
             5.0050E-01
 GRADIENT:  -6.0131E-02  2.0213E+01  5.6258E+00 -1.2786E+02  1.5640E+02  1.7012E-01  2.6638E+01  0.0000E+00 -1.8473E+00 -9.8158E-01
             1.6512E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1719
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.5679E-04  1.5308E-02 -1.1919E-04 -8.4049E-03 -9.5079E-03
 SE:             2.9656E-02  9.5112E-03  2.0125E-04  2.8683E-02  2.3751E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7695E-01  1.0752E-01  5.5368E-01  7.6950E-01  6.8892E-01

 ETASHRINKSD(%)  6.4740E-01  6.8136E+01  9.9326E+01  3.9077E+00  2.0432E+01
 ETASHRINKVR(%)  1.2906E+00  8.9847E+01  9.9995E+01  7.6626E+00  3.6689E+01
 EBVSHRINKSD(%)  7.4068E-01  8.0208E+01  9.9242E+01  3.6677E+00  1.8973E+01
 EBVSHRINKVR(%)  1.4759E+00  9.6083E+01  9.9994E+01  7.2008E+00  3.4346E+01
 RELATIVEINF(%)  9.8261E+01  2.0487E+00  5.4271E-04  4.9427E+01  6.2497E+00
 EPSSHRINKSD(%)  3.0059E+01
 EPSSHRINKVR(%)  5.1083E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1899.8485225052818     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -980.90998930060914     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1899.849       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.16E-01  3.92E-02  8.57E-01  1.55E+00  5.81E-01  1.02E+00  6.27E+00  1.00E-02  8.25E-01  7.93E-01  1.49E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,      143.775
Stop Time:
Sat Oct 23 22:01:11 CDT 2021
