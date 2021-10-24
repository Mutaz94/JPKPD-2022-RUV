Sun Oct 24 01:23:48 CDT 2021
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
$DATA ../../../../data/SD3/D2/dat89.csv ignore=@
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
Current Date:       24 OCT 2021
Days until program expires : 175
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2019.05959590888        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.6468E+02 -1.1830E+02 -7.8979E+01 -1.3605E+02  1.3898E+02 -1.4010E+01 -6.1400E+01  5.7819E+00 -1.6398E+02 -1.0849E+01
             5.6614E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2079.68598632631        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      130
 NPARAMETR:  1.0669E+00  1.0002E+00  1.0897E+00  1.0714E+00  9.3420E-01  1.1137E+00  1.1073E+00  9.8524E-01  1.3388E+00  9.9046E-01
             9.1687E-01
 PARAMETER:  1.6476E-01  1.0020E-01  1.8589E-01  1.6901E-01  3.1938E-02  2.0768E-01  2.0195E-01  8.5132E-02  3.9180E-01  9.0415E-02
             1.3213E-02
 GRADIENT:   1.4491E+01 -3.1644E+01 -2.0412E+01 -6.2921E+01  7.7301E-01 -7.7329E+00 -1.7767E+01  6.8532E-01 -3.3601E+01 -1.8914E+00
            -4.5160E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2086.22452454890        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0619E+00  6.1995E-01  1.9331E+00  1.4507E+00  1.0846E+00  1.0753E+00  1.6840E+00  1.3773E+00  1.2232E+00  1.1848E+00
             9.2363E-01
 PARAMETER:  1.6004E-01 -3.7812E-01  7.5910E-01  4.7205E-01  1.8121E-01  1.7265E-01  6.2116E-01  4.2014E-01  3.0150E-01  2.6957E-01
             2.0553E-02
 GRADIENT:   1.5872E+01  1.0641E+01  1.4430E+01  2.4143E+01  1.8696E+01 -2.0872E+01  8.3869E-01 -1.7254E+01 -7.1267E+00  3.6688E+00
             2.8764E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2089.43442834286        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      484
 NPARAMETR:  1.0514E+00  6.3483E-01  1.7803E+00  1.4097E+00  1.0204E+00  1.1244E+00  1.4801E+00  1.4982E+00  1.2805E+00  1.0939E+00
             9.1783E-01
 PARAMETER:  1.5014E-01 -3.5439E-01  6.7679E-01  4.4338E-01  1.2024E-01  2.1728E-01  4.9214E-01  5.0427E-01  3.4725E-01  1.8975E-01
             1.4255E-02
 GRADIENT:  -3.4453E+00  7.0344E+00  6.7245E+00  2.6154E+00 -1.0273E+01 -1.7981E+00  1.6008E+00 -7.7784E-01  1.6622E+00  1.0262E+00
             3.4545E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2090.72367754427        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      661
 NPARAMETR:  1.0520E+00  4.0637E-01  1.7871E+00  1.5689E+00  9.7900E-01  1.1387E+00  7.9542E-01  1.5509E+00  1.2198E+00  1.0940E+00
             9.1482E-01
 PARAMETER:  1.5073E-01 -8.0049E-01  6.8058E-01  5.5040E-01  7.8777E-02  2.2992E-01 -1.2888E-01  5.3882E-01  2.9866E-01  1.8985E-01
             1.0967E-02
 GRADIENT:   2.8035E+00  4.5915E+00 -7.8840E+00  2.0039E+01  2.8802E+00  4.5099E+00  6.0783E-01  1.9590E+00  3.8335E+00  3.7425E+00
            -2.1728E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2091.51121129017        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      838
 NPARAMETR:  1.0450E+00  2.2750E-01  1.9604E+00  1.6746E+00  9.6661E-01  1.1262E+00  3.0641E-01  1.6655E+00  1.1188E+00  1.0562E+00
             9.1827E-01
 PARAMETER:  1.4400E-01 -1.3806E+00  7.7313E-01  6.1558E-01  6.6036E-02  2.1887E-01 -1.0828E+00  6.1014E-01  2.1227E-01  1.5467E-01
             1.4735E-02
 GRADIENT:  -5.1500E+00  2.1644E+00 -2.8196E-02  2.6193E+00 -3.1928E-01  6.9546E-01  4.3040E-02 -2.5078E-01 -1.4834E+00 -6.3477E-01
             5.6365E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2091.51277679469        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1017
 NPARAMETR:  1.0449E+00  2.0582E-01  1.9814E+00  1.6903E+00  9.6543E-01  1.1254E+00  2.5163E-01  1.6822E+00  1.1087E+00  1.0552E+00
             9.1823E-01
 PARAMETER:  1.4397E-01 -1.4807E+00  7.8380E-01  6.2492E-01  6.4816E-02  2.1816E-01 -1.2798E+00  6.2012E-01  2.0322E-01  1.5376E-01
             1.4690E-02
 GRADIENT:  -4.6872E+00  2.0762E+00  3.4608E-01  2.8788E+00 -8.5493E-01  4.6858E-01  2.4695E-02 -2.4868E-01 -1.4939E+00 -5.8927E-01
             5.4800E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2091.51353798348        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1197
 NPARAMETR:  1.0449E+00  1.9496E-01  1.9901E+00  1.6980E+00  9.6464E-01  1.1252E+00  2.2549E-01  1.6895E+00  1.1038E+00  1.0547E+00
             9.1819E-01
 PARAMETER:  1.4396E-01 -1.5350E+00  7.8817E-01  6.2945E-01  6.4001E-02  2.1795E-01 -1.3895E+00  6.2442E-01  1.9880E-01  1.5324E-01
             1.4648E-02
 GRADIENT:  -4.4232E+00  1.9952E+00  3.8659E-01  2.8462E+00 -9.1313E-01  4.1623E-01  1.8165E-02 -2.3779E-01 -1.4583E+00 -5.5951E-01
             5.2395E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2091.51368688261        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1377
 NPARAMETR:  1.0449E+00  1.8830E-01  1.9952E+00  1.7027E+00  9.6415E-01  1.1251E+00  2.1006E-01  1.6938E+00  1.1009E+00  1.0543E+00
             9.1816E-01
 PARAMETER:  1.4395E-01 -1.5697E+00  7.9073E-01  6.3219E-01  6.3487E-02  2.1784E-01 -1.4604E+00  6.2697E-01  1.9612E-01  1.5292E-01
             1.4620E-02
 GRADIENT:  -4.2645E+00  1.9394E+00  3.9179E-01  2.7951E+00 -9.1751E-01  3.9159E-01  1.4902E-02 -2.3040E-01 -1.4249E+00 -5.4110E-01
             5.0762E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2091.65594245828        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1566
 NPARAMETR:  1.0494E+00  1.7200E-01  1.9967E+00  1.6932E+00  9.6097E-01  1.1256E+00  3.0826E-02  1.6970E+00  1.0959E+00  1.0576E+00
             9.1730E-01
 PARAMETER:  1.4821E-01 -1.6603E+00  7.9150E-01  6.2660E-01  6.0186E-02  2.1829E-01 -3.3794E+00  6.2885E-01  1.9157E-01  1.5600E-01
             1.3675E-02
 GRADIENT:   3.9556E+00  2.4387E-01  5.7346E-01 -1.8643E+01 -4.5173E-01  6.4992E-01  3.9570E-04  4.6028E-02 -3.7179E-01  2.6012E-01
             1.5406E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2091.65808611390        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1751             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0492E+00  1.7176E-01  1.9910E+00  1.6939E+00  9.6054E-01  1.1255E+00  1.0000E-02  1.6952E+00  1.0976E+00  1.0554E+00
             9.1732E-01
 PARAMETER:  1.4802E-01 -1.6617E+00  7.8863E-01  6.2705E-01  5.9740E-02  2.1826E-01 -5.6800E+00  6.2778E-01  1.9309E-01  1.5389E-01
             1.3702E-02
 GRADIENT:   7.0073E+02  2.2073E+01  1.4976E+01  1.2416E+03  8.3816E+00  1.1609E+02  0.0000E+00  5.0902E+00  4.4512E+01  1.7200E+00
             9.9635E-01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -2091.65808611390        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     1816
 NPARAMETR:  1.0492E+00  1.7245E-01  1.9907E+00  1.6943E+00  9.5958E-01  1.1255E+00  1.0000E-02  1.6934E+00  1.0972E+00  1.0555E+00
             9.1733E-01
 PARAMETER:  1.4802E-01 -1.6617E+00  7.8863E-01  6.2705E-01  5.9740E-02  2.1826E-01 -5.6800E+00  6.2778E-01  1.9309E-01  1.5389E-01
             1.3702E-02
 GRADIENT:   2.3232E-02 -3.0695E-02  3.2816E-02 -2.6041E-01  3.5299E-01  1.8335E-03  0.0000E+00  8.9512E-02  6.9812E-02 -1.6555E-02
            -3.0173E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1816
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.3440E-05 -5.7361E-05 -3.9831E-02 -4.3577E-03 -4.4404E-02
 SE:             2.9922E-02  2.1491E-05  1.9425E-02  2.9759E-02  2.0824E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9884E-01  7.6052E-03  4.0314E-02  8.8358E-01  3.2984E-02

 ETASHRINKSD(%)  1.0000E-10  9.9928E+01  3.4925E+01  3.0386E-01  3.0235E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  5.7652E+01  6.0679E-01  5.1329E+01
 EBVSHRINKSD(%)  2.1623E-01  9.9929E+01  3.7758E+01  7.8763E-01  2.6860E+01
 EBVSHRINKVR(%)  4.3198E-01  1.0000E+02  6.1259E+01  1.5691E+00  4.6505E+01
 RELATIVEINF(%)  9.8495E+01  6.7073E-06  1.4646E+01  1.5912E+01  1.5895E+01
 EPSSHRINKSD(%)  3.5294E+01
 EPSSHRINKVR(%)  5.8131E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2091.6580861138982     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1172.7195529092255     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2091.658       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.72E-01  1.99E+00  1.69E+00  9.61E-01  1.13E+00  1.00E-02  1.70E+00  1.10E+00  1.06E+00  9.17E-01
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       61.056
Stop Time:
Sun Oct 24 01:24:00 CDT 2021
