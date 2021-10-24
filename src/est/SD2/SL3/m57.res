Sat Oct 23 19:24:17 CDT 2021
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
$DATA ../../../../data/SD2/SL3/dat57.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      798
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

 TOT. NO. OF OBS RECS:      698
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2156.14943882307        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6054E+02 -1.0327E+02 -2.3268E+01  9.3647E+01  1.4246E+02  7.0491E+01 -5.2543E+01 -1.1705E+01  6.7429E+00 -5.1977E+00
            -1.6566E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2580.59769382676        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0693E+00  1.2425E+00  9.9218E-01  9.0410E-01  1.0461E+00  8.9932E-01  1.1008E+00  8.8305E-01  7.9209E-01  9.9368E-01
             1.8161E+00
 PARAMETER:  1.6703E-01  3.1716E-01  9.2154E-02 -8.1274E-04  1.4509E-01 -6.1218E-03  1.9602E-01 -2.4369E-02 -1.3307E-01  9.3656E-02
             6.9667E-01
 GRADIENT:   2.8682E+02  4.6604E+01 -1.9309E+00  1.3000E+00  1.2466E+01  9.6540E+00  1.2846E+01  3.3333E+00 -1.3022E+01 -1.4532E+01
             1.4924E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2585.53709273647        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      266
 NPARAMETR:  1.0651E+00  1.5679E+00  1.0599E+00  7.6800E-01  1.2942E+00  9.4405E-01  8.2468E-01  7.7371E-01  1.0848E+00  1.2698E+00
             1.8034E+00
 PARAMETER:  1.6302E-01  5.4972E-01  1.5820E-01 -1.6397E-01  3.5791E-01  4.2420E-02 -9.2766E-02 -1.5656E-01  1.8137E-01  3.3884E-01
             6.8966E-01
 GRADIENT:   4.2325E+01  4.6438E+01 -7.4767E+00  6.9263E+01  3.0889E+01  1.7317E+01 -1.0142E+00 -7.6123E-01 -1.9506E+00  6.7709E-01
             2.2292E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2589.79808881750        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.0347E+00  1.6635E+00  8.6900E-01  6.6595E-01  1.2520E+00  8.8267E-01  7.6945E-01  5.3626E-01  1.2051E+00  1.2903E+00
             1.7990E+00
 PARAMETER:  1.3413E-01  6.0894E-01 -4.0415E-02 -3.0654E-01  3.2475E-01 -2.4809E-02 -1.6207E-01 -5.2313E-01  2.8658E-01  3.5485E-01
             6.8721E-01
 GRADIENT:  -3.2484E+01  1.5597E+01  1.8187E+00  1.6227E+01 -1.3829E+01 -7.2623E+00  1.1145E+00  4.7646E-02  7.4930E-01  5.8634E+00
             4.1888E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2594.37539083248        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  1.0511E+00  2.1955E+00  5.2245E-01  3.2436E-01  1.5729E+00  8.9259E-01  6.3392E-01  2.5556E-01  1.9982E+00  1.4861E+00
             1.7694E+00
 PARAMETER:  1.4984E-01  8.8639E-01 -5.4923E-01 -1.0259E+00  5.5291E-01 -1.3627E-02 -3.5583E-01 -1.2643E+00  7.9225E-01  4.9616E-01
             6.7064E-01
 GRADIENT:   1.1705E+01  3.9548E+01 -2.0900E+00  1.7461E+01  7.6236E+00 -3.1306E+00 -1.8251E+00  5.8663E-02  8.2390E-01  2.5443E+00
            -1.6341E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2595.66064005380        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      796
 NPARAMETR:  1.0463E+00  2.3755E+00  3.7393E-01  1.8463E-01  1.6889E+00  9.0047E-01  6.2145E-01  9.3680E-02  2.6550E+00  1.5555E+00
             1.7872E+00
 PARAMETER:  1.4523E-01  9.6522E-01 -8.8369E-01 -1.5894E+00  6.2408E-01 -4.8338E-03 -3.7570E-01 -2.2679E+00  1.0765E+00  5.4179E-01
             6.8063E-01
 GRADIENT:  -3.5827E-01 -1.4753E+00 -1.8939E+00  2.4263E+00  3.2868E+00  6.7856E-01  1.5463E-01  1.4055E-02  3.5769E-01  8.5324E-02
             1.3483E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2595.73899732006        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      976
 NPARAMETR:  1.0463E+00  2.3776E+00  3.9556E-01  1.8046E-01  1.6882E+00  8.9883E-01  6.1862E-01  3.1543E-02  2.7023E+00  1.5577E+00
             1.7852E+00
 PARAMETER:  1.4524E-01  9.6610E-01 -8.2745E-01 -1.6123E+00  6.2369E-01 -6.6661E-03 -3.8027E-01 -3.3564E+00  1.0941E+00  5.4324E-01
             6.7953E-01
 GRADIENT:  -1.4043E-01 -3.0694E+00  2.4919E-01  2.6873E-01 -6.5613E-01 -1.1175E-02  5.0868E-02  1.3096E-03 -7.5452E-01 -2.6563E-01
            -2.6539E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2595.74996066877        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1160
 NPARAMETR:  1.0468E+00  2.3736E+00  3.9502E-01  1.8037E-01  1.6889E+00  8.9887E-01  6.1843E-01  1.0000E-02  2.7161E+00  1.5595E+00
             1.7853E+00
 PARAMETER:  1.4570E-01  9.6440E-01 -8.2883E-01 -1.6127E+00  6.2408E-01 -6.6152E-03 -3.8057E-01 -5.5248E+00  1.0992E+00  5.4434E-01
             6.7960E-01
 GRADIENT:   1.2237E+00 -1.2072E+01 -1.5340E-01  8.2343E-02 -2.8005E-02  2.9531E-02  1.0664E-01  0.0000E+00  2.9696E-01  1.3812E-01
             1.2619E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -2595.74996066877        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1182
 NPARAMETR:  1.0468E+00  2.3736E+00  3.9502E-01  1.8037E-01  1.6889E+00  8.9887E-01  6.1843E-01  1.0000E-02  2.7161E+00  1.5595E+00
             1.7853E+00
 PARAMETER:  1.4570E-01  9.6440E-01 -8.2883E-01 -1.6127E+00  6.2408E-01 -6.6152E-03 -3.8057E-01 -5.5248E+00  1.0992E+00  5.4434E-01
             6.7960E-01
 GRADIENT:   1.2237E+00 -1.2072E+01 -1.5340E-01  8.2343E-02 -2.8005E-02  2.9531E-02  1.0664E-01  0.0000E+00  2.9696E-01  1.3812E-01
             1.2619E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1182
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.1697E-04 -3.9319E-02 -6.6218E-05  4.7394E-02 -3.2786E-02
 SE:             2.9652E-02  2.4394E-02  3.0687E-05  1.9017E-02  2.5354E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8071E-01  1.0699E-01  3.0938E-02  1.2694E-02  1.9596E-01

 ETASHRINKSD(%)  6.6065E-01  1.8277E+01  9.9897E+01  3.6292E+01  1.5061E+01
 ETASHRINKVR(%)  1.3169E+00  3.3214E+01  1.0000E+02  5.9413E+01  2.7853E+01
 EBVSHRINKSD(%)  9.9585E-01  1.4772E+01  9.9887E+01  4.5919E+01  1.0790E+01
 EBVSHRINKVR(%)  1.9818E+00  2.7361E+01  1.0000E+02  7.0753E+01  2.0416E+01
 RELATIVEINF(%)  9.7982E+01  1.6072E+01  7.2819E-05  6.1117E+00  5.1543E+01
 EPSSHRINKSD(%)  2.2212E+01
 EPSSHRINKVR(%)  3.9490E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          698
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1282.8381923537231     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2595.7499606687747     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1312.9117683150516     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2595.750       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.37E+00  3.95E-01  1.80E-01  1.69E+00  8.99E-01  6.18E-01  1.00E-02  2.72E+00  1.56E+00  1.79E+00
 


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
 #CPUT: Total CPU Time in Seconds,      106.579
Stop Time:
Sat Oct 23 19:24:33 CDT 2021
