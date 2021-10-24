Sat Oct 23 15:11:18 CDT 2021
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
$DATA ../../../../data/SD1/SL2/dat61.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2737.84600416361        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0099E+02  1.0722E+02  9.0666E+01  7.3243E+01  6.7156E+01  6.3548E+01 -6.7610E+01 -1.2757E+02 -2.8479E+01  1.0702E+01
            -2.0804E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3242.06581844550        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9547E-01  9.6108E-01  9.3393E-01  1.0288E+00  9.4363E-01  9.1534E-01  1.1665E+00  1.1539E+00  9.1030E-01  8.8328E-01
             1.7670E+00
 PARAMETER:  9.5461E-02  6.0304E-02  3.1647E-02  1.2840E-01  4.1980E-02  1.1542E-02  2.5399E-01  2.4315E-01  6.0228E-03 -2.4108E-02
             6.6929E-01
 GRADIENT:   5.9318E+01  7.9268E+00  5.2585E+00  1.1389E+01 -3.2458E+00 -6.3745E+00  3.0512E-01 -1.0524E+00 -7.1014E+00  2.9485E+00
             2.1521E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3245.27524425943        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      263
 NPARAMETR:  1.0041E+00  1.1762E+00  1.2041E+00  9.3312E-01  1.2090E+00  9.4345E-01  1.0717E+00  1.8936E+00  1.0029E+00  9.7544E-01
             1.7971E+00
 PARAMETER:  1.0409E-01  2.6230E-01  2.8572E-01  3.0778E-02  2.8981E-01  4.1790E-02  1.6923E-01  7.3850E-01  1.0287E-01  7.5134E-02
             6.8615E-01
 GRADIENT:  -6.0216E+01 -1.7896E+01 -8.5291E+00 -5.6230E-01  1.3733E+01 -1.0087E+01  3.3459E+00  5.1733E+00  1.3636E+01 -9.0458E+00
             5.3538E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3248.36098312108        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.0352E+00  1.2327E+00  1.2850E+00  9.1077E-01  1.2663E+00  9.6362E-01  1.0196E+00  1.9947E+00  9.0386E-01  1.1007E+00
             1.7605E+00
 PARAMETER:  1.3456E-01  3.0920E-01  3.5075E-01  6.5302E-03  3.3613E-01  6.2941E-02  1.1942E-01  7.9048E-01 -1.0856E-03  1.9590E-01
             6.6559E-01
 GRADIENT:   1.5001E+01 -7.9493E+00 -2.2078E+00  5.2096E+00  6.7184E-01 -2.5806E-01 -2.2586E+00 -2.8029E-01 -1.1105E-01  2.2749E+00
             2.5951E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3249.58228749529        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      619
 NPARAMETR:  1.0249E+00  1.3465E+00  1.4600E+00  8.5200E-01  1.4025E+00  9.7585E-01  1.0190E+00  2.5688E+00  8.5707E-01  1.1519E+00
             1.7529E+00
 PARAMETER:  1.2458E-01  3.9748E-01  4.7844E-01 -6.0168E-02  4.3823E-01  7.5554E-02  1.1880E-01  1.0434E+00 -5.4236E-02  2.4143E-01
             6.6126E-01
 GRADIENT:  -8.4885E+00  6.7036E+00 -4.2795E+00  5.7398E+00  5.0283E+00  4.3849E+00 -2.0846E-01 -6.8686E-01  3.0111E+00 -1.4725E+00
             1.8198E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3249.69472623553        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      801
 NPARAMETR:  1.0315E+00  1.3488E+00  1.4685E+00  8.5134E-01  1.4064E+00  9.6883E-01  1.0383E+00  2.5690E+00  8.0721E-01  1.1511E+00
             1.7571E+00
 PARAMETER:  1.3106E-01  3.9919E-01  4.8426E-01 -6.0948E-02  4.4101E-01  6.8329E-02  1.3759E-01  1.0435E+00 -1.1418E-01  2.4076E-01
             6.6368E-01
 GRADIENT:   6.3961E+00  9.0591E+00 -3.7983E+00  5.9927E+00  4.4714E+00  1.7306E+00 -4.1752E-01 -2.2712E+00 -2.3724E-01 -2.9430E+00
             5.0780E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3249.76687873108        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      981             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0287E+00  1.3467E+00  1.4799E+00  8.5078E-01  1.4057E+00  9.6449E-01  1.0427E+00  2.5824E+00  8.0859E-01  1.1510E+00
             1.7512E+00
 PARAMETER:  1.2830E-01  3.9767E-01  4.9201E-01 -6.1601E-02  4.4053E-01  6.3840E-02  1.4178E-01  1.0487E+00 -1.1246E-01  2.4061E-01
             6.6031E-01
 GRADIENT:   1.9233E+02  1.4603E+02  3.1805E+00  2.5788E+01  5.9058E+01  1.6631E+01  4.7795E+00  3.9067E+00  1.6018E+00  1.0496E+00
             1.0671E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3249.84623712997        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1158
 NPARAMETR:  1.0284E+00  1.3344E+00  1.5080E+00  8.5714E-01  1.4041E+00  9.6411E-01  1.0469E+00  2.5729E+00  8.0093E-01  1.1627E+00
             1.7553E+00
 PARAMETER:  1.2796E-01  3.8850E-01  5.1078E-01 -5.4155E-02  4.3941E-01  6.3454E-02  1.4587E-01  1.0450E+00 -1.2199E-01  2.5074E-01
             6.6265E-01
 GRADIENT:  -7.6034E-01  3.7500E+00 -1.0314E+00  1.1659E+00  3.6858E+00 -1.1486E-01 -2.8346E-01 -3.7986E+00 -4.0326E-03 -7.6228E-02
             4.8850E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3250.06300906723        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1328
 NPARAMETR:  1.0280E+00  1.2704E+00  1.5831E+00  8.9963E-01  1.3720E+00  9.6406E-01  1.0915E+00  2.6041E+00  7.7755E-01  1.1247E+00
             1.7446E+00
 PARAMETER:  1.2744E-01  3.3916E-01  5.5912E-01 -5.7239E-03  4.1613E-01  6.3322E-02  1.8803E-01  1.0575E+00 -1.5314E-01  2.1760E-01
             6.5625E-01
 GRADIENT:  -1.8340E+00 -1.8571E+02 -1.1589E+02  3.2762E+02 -1.4840E+02 -1.8433E-01  5.1290E-01  5.8640E+01 -6.8312E-01  2.9733E+02
            -1.0371E+02
 NUMSIGDIG:         1.9         2.3         2.3         2.3         2.3         2.0         1.5         2.3         0.9         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1328
 NO. OF SIG. DIGITS IN FINAL EST.:  0.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0761E-03 -1.7389E-02 -3.1660E-02  1.1688E-02 -3.7748E-02
 SE:             2.9756E-02  2.4100E-02  2.1744E-02  2.1418E-02  2.2697E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4437E-01  4.7058E-01  1.4538E-01  5.8525E-01  9.6291E-02

 ETASHRINKSD(%)  3.1524E-01  1.9263E+01  2.7156E+01  2.8248E+01  2.3962E+01
 ETASHRINKVR(%)  6.2949E-01  3.4815E+01  4.6937E+01  4.8517E+01  4.2182E+01
 EBVSHRINKSD(%)  7.6244E-01  1.9462E+01  3.1480E+01  3.1874E+01  2.1442E+01
 EBVSHRINKVR(%)  1.5191E+00  3.5136E+01  5.3050E+01  5.3588E+01  3.8286E+01
 RELATIVEINF(%)  9.8466E+01  1.6857E+01  3.3425E+01  1.1881E+01  2.7321E+01
 EPSSHRINKSD(%)  1.9606E+01
 EPSSHRINKVR(%)  3.5368E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3250.0630090672325     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1595.9736492988218     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3250.063       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.27E+00  1.58E+00  9.00E-01  1.37E+00  9.64E-01  1.09E+00  2.61E+00  7.76E-01  1.12E+00  1.74E+00
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,      103.665
Stop Time:
Sat Oct 23 15:11:35 CDT 2021
