Sun Oct 24 04:02:56 CDT 2021
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
$DATA ../../../../data/SD4/TD2/dat41.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1702.02332504103        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5053E+02 -4.5389E+01  1.5116E-01 -3.1176E+01  3.8784E+01  4.8093E+01  1.5189E+00  3.1020E-01  2.8370E+01 -1.4440E+01
             2.4496E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1711.24456721683        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0364E+00  1.0637E+00  9.3263E-01  1.0396E+00  9.6985E-01  9.6703E-01  1.0074E+00  9.9829E-01  8.2522E-01  1.1035E+00
             9.2511E-01
 PARAMETER:  1.3575E-01  1.6176E-01  3.0257E-02  1.3886E-01  6.9385E-02  6.6476E-02  1.0741E-01  9.8291E-02 -9.2103E-02  1.9853E-01
             2.2160E-02
 GRADIENT:   5.2671E+00  8.6019E+00 -3.4339E+00  1.6708E+01 -3.2681E-01 -3.0631E+00 -6.3650E+00  4.3190E+00 -4.6832E+00  4.1264E+00
            -6.2233E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1712.43796807161        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0344E+00  8.9140E-01  8.0450E-01  1.1411E+00  8.1987E-01  9.5977E-01  1.3593E+00  6.3734E-01  7.1930E-01  1.0132E+00
             9.2495E-01
 PARAMETER:  1.3383E-01 -1.4966E-02 -1.1753E-01  2.3198E-01 -9.8606E-02  5.8936E-02  4.0700E-01 -3.5046E-01 -2.2948E-01  1.1307E-01
             2.1982E-02
 GRADIENT:  -1.0445E+00  2.1651E+01 -1.2413E+01  4.0202E+01  9.4148E-01 -6.4976E+00  4.4118E+00  3.8554E+00 -4.5957E+00  7.2521E+00
            -2.6334E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1714.09040740295        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0356E+00  8.0559E-01  7.2100E-01  1.1527E+00  7.4093E-01  9.8405E-01  1.3906E+00  2.9740E-01  7.3469E-01  9.3961E-01
             9.3030E-01
 PARAMETER:  1.3501E-01 -1.1617E-01 -2.2712E-01  2.4211E-01 -1.9985E-01  8.3921E-02  4.2975E-01 -1.1127E+00 -2.0831E-01  3.7706E-02
             2.7749E-02
 GRADIENT:   7.1847E-01 -5.5459E+00 -4.5892E+00 -1.1810E+01  2.9714E+00  3.2250E+00  9.5937E-01  7.0592E-01  1.1095E+00  1.3831E+00
             1.1217E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1714.22170156542        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.0347E+00  7.2998E-01  7.3993E-01  1.2029E+00  7.2385E-01  9.7376E-01  1.4892E+00  2.1088E-01  7.1511E-01  9.6144E-01
             9.2866E-01
 PARAMETER:  1.3409E-01 -2.1474E-01 -2.0120E-01  2.8477E-01 -2.2317E-01  7.3414E-02  4.9827E-01 -1.4565E+00 -2.3532E-01  6.0674E-02
             2.5983E-02
 GRADIENT:   1.8898E-01 -9.4882E-01 -1.0541E+00 -9.0612E-01  6.7534E-01 -6.6923E-01  6.6028E-02  1.7408E-01  2.3833E-02  4.8801E-01
             2.6136E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1714.23721136917        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0346E+00  7.1444E-01  7.3251E-01  1.2094E+00  7.1414E-01  9.7679E-01  1.5104E+00  1.3869E-01  7.1280E-01  9.5903E-01
             9.2812E-01
 PARAMETER:  1.3403E-01 -2.3626E-01 -2.1128E-01  2.9015E-01 -2.3668E-01  7.6517E-02  5.1238E-01 -1.8755E+00 -2.3856E-01  5.8165E-02
             2.5404E-02
 GRADIENT:   1.1949E-01 -1.9528E+00 -4.0946E-01 -2.7905E+00  6.7976E-01  5.4858E-01 -1.9321E-02  6.0398E-02  2.2738E-01  8.0848E-02
             1.6978E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1714.36151591284        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  1.0353E+00  7.8226E-01  7.1353E-01  1.1705E+00  7.2647E-01  9.7614E-01  1.4131E+00  1.3404E-02  7.2725E-01  9.5066E-01
             9.2850E-01
 PARAMETER:  1.3471E-01 -1.4557E-01 -2.3753E-01  2.5742E-01 -2.1956E-01  7.5848E-02  4.4582E-01 -4.2122E+00 -2.1848E-01  4.9401E-02
             2.5820E-02
 GRADIENT:  -8.6913E-03  7.7167E-02  1.7913E-02  6.7650E-02 -8.1323E-02  1.9250E-02  3.2591E-03  8.1900E-04 -1.0002E-02  1.5842E-02
            -8.5856E-03

0ITERATION NO.:   33    OBJECTIVE VALUE:  -1714.36406167293        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:     1170
 NPARAMETR:  1.0367E+00  7.8216E-01  7.1351E-01  1.1698E+00  7.2647E-01  9.7631E-01  1.4142E+00  1.0000E-02  7.2743E-01  9.5059E-01
             9.2852E-01
 PARAMETER:  1.3601E-01 -1.4570E-01 -2.3756E-01  2.5683E-01 -2.1955E-01  7.6024E-02  4.4654E-01 -4.6909E+00 -2.1824E-01  4.9324E-02
             2.5839E-02
 GRADIENT:   3.0693E+00 -4.2245E-01  3.7710E-01 -1.6870E+00 -3.6151E-01  9.5733E-02  9.6447E-02  0.0000E+00  6.8349E-02  5.1853E-02
            -5.1857E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1170
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4478E-05  7.5245E-03 -5.3613E-04 -9.6938E-03 -5.2911E-03
 SE:             2.9854E-02  2.0675E-02  2.2005E-04  2.4369E-02  2.4374E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9961E-01  7.1590E-01  1.4834E-02  6.9079E-01  8.2815E-01

 ETASHRINKSD(%)  1.0000E-10  3.0737E+01  9.9263E+01  1.8360E+01  1.8344E+01
 ETASHRINKVR(%)  1.0000E-10  5.2027E+01  9.9995E+01  3.3348E+01  3.3323E+01
 EBVSHRINKSD(%)  3.8343E-01  3.1216E+01  9.9349E+01  1.8358E+01  1.6479E+01
 EBVSHRINKVR(%)  7.6539E-01  5.2688E+01  9.9996E+01  3.3346E+01  3.0242E+01
 RELATIVEINF(%)  9.8778E+01  3.5597E+00  5.0587E-04  5.8895E+00  5.8184E+00
 EPSSHRINKSD(%)  4.4613E+01
 EPSSHRINKVR(%)  6.9323E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1714.3640616729263     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -979.21323510918808     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1714.364       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  7.82E-01  7.14E-01  1.17E+00  7.26E-01  9.76E-01  1.41E+00  1.00E-02  7.27E-01  9.51E-01  9.29E-01
 


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
 #CPUT: Total CPU Time in Seconds,       32.702
Stop Time:
Sun Oct 24 04:03:04 CDT 2021
