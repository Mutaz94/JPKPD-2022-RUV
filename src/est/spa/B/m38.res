Wed Sep 29 11:11:37 CDT 2021
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
$DATA ../../../../data/spa/B/dat38.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1618.54229434050        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6499E+02  2.0800E+01 -1.5980E+01  5.0115E+01 -1.9552E+01  3.0379E+01 -1.9661E+01  1.3128E+01 -3.7576E+01  1.1045E+01
            -3.4690E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1629.88456009477        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0010E+00  1.0519E+00  1.1440E+00  9.9275E-01  1.0806E+00  1.0899E+00  1.1144E+00  9.3147E-01  1.2002E+00  9.6022E-01
             1.1061E+00
 PARAMETER:  1.0099E-01  1.5064E-01  2.3452E-01  9.2723E-02  1.7749E-01  1.8609E-01  2.0829E-01  2.9009E-02  2.8245E-01  5.9403E-02
             2.0080E-01
 GRADIENT:   1.5458E+01  3.2205E+00  4.0364E+00  2.2100E+00 -2.7392E+01  8.0172E+00 -3.2230E+00  4.8757E+00  1.9110E+00 -1.8480E+00
             7.7754E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1631.25771011880        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  1.0035E+00  1.1329E+00  1.1340E+00  9.4336E-01  1.1520E+00  1.0702E+00  1.1072E+00  7.9052E-01  1.2646E+00  1.0781E+00
             1.0686E+00
 PARAMETER:  1.0349E-01  2.2479E-01  2.2571E-01  4.1696E-02  2.4147E-01  1.6787E-01  2.0187E-01 -1.3506E-01  3.3478E-01  1.7523E-01
             1.6637E-01
 GRADIENT:   2.0740E+01 -2.4288E+00 -3.0166E+00  7.5539E+00  6.3764E-01  4.7799E-01  1.5113E+00  2.2691E+00  4.8296E+00  4.4566E+00
            -4.7661E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1633.24017146153        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.8648E-01  1.3666E+00  8.9163E-01  7.7618E-01  1.1692E+00  1.0701E+00  9.9430E-01  2.8964E-01  1.3948E+00  1.0131E+00
             1.0951E+00
 PARAMETER:  8.6392E-02  4.1231E-01 -1.4705E-02 -1.5337E-01  2.5629E-01  1.6776E-01  9.4286E-02 -1.1391E+00  4.3278E-01  1.1301E-01
             1.9088E-01
 GRADIENT:  -1.6292E+01 -9.7627E+00 -3.5966E+00  1.1856E+00  4.4196E+00 -8.2216E-02 -1.7012E+00  3.7590E-01 -4.9868E-01 -1.0238E-01
             4.7635E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1634.16244082073        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  1.0040E+00  1.6959E+00  7.5403E-01  5.6941E-01  1.2985E+00  1.0708E+00  8.6379E-01  7.7875E-02  1.7638E+00  1.0878E+00
             1.0799E+00
 PARAMETER:  1.0399E-01  6.2819E-01 -1.8233E-01 -4.6315E-01  3.6122E-01  1.6843E-01 -4.6423E-02 -2.4526E+00  6.6747E-01  1.8412E-01
             1.7685E-01
 GRADIENT:   1.6228E+01  4.9444E+00  2.0816E-01  5.8293E+00 -2.1125E+00 -8.4836E-02  7.7046E-02  2.6088E-02  1.4693E-01  3.4839E-01
            -4.0037E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1634.16602838871        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  1.0033E+00  1.7365E+00  7.3642E-01  5.4397E-01  1.3160E+00  1.0706E+00  8.4948E-01  6.1050E-02  1.8208E+00  1.0973E+00
             1.0820E+00
 PARAMETER:  1.0325E-01  6.5185E-01 -2.0596E-01 -5.0887E-01  3.7463E-01  1.6824E-01 -6.3136E-02 -2.6961E+00  6.9930E-01  1.9286E-01
             1.7881E-01
 GRADIENT:   1.4587E+01  7.1994E+00  3.8795E-01  6.2920E+00 -2.2421E+00 -1.2709E-01 -1.5004E-01  1.5720E-02 -5.9593E-02  2.6893E-01
            -3.7067E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1634.16710619579        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0028E+00  1.7536E+00  7.2850E-01  5.3260E-01  1.3237E+00  1.0706E+00  8.4389E-01  5.4233E-02  1.8471E+00  1.1013E+00
             1.0831E+00
 PARAMETER:  1.0281E-01  6.6166E-01 -2.1677E-01 -5.2998E-01  3.8042E-01  1.6818E-01 -6.9737E-02 -2.8145E+00  7.1361E-01  1.9651E-01
             1.7979E-01
 GRADIENT:   1.3667E+01  7.3103E+00  4.0776E-01  6.1257E+00 -2.1709E+00 -1.2998E-01 -1.8361E-01  1.2351E-02 -9.4900E-02  2.4087E-01
            -3.4963E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1634.30562666163        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1248
 NPARAMETR:  9.9658E-01  1.7443E+00  7.2811E-01  5.2614E-01  1.3257E+00  1.0721E+00  8.4230E-01  1.0941E-02  1.8558E+00  1.0957E+00
             1.0911E+00
 PARAMETER:  9.6570E-02  6.5634E-01 -2.1730E-01 -5.4219E-01  3.8195E-01  1.6958E-01 -7.1617E-02 -4.4153E+00  7.1832E-01  1.9140E-01
             1.8717E-01
 GRADIENT:   1.6316E+00 -1.0506E+01 -6.9349E-01 -3.3224E-01  1.2009E+00  5.9742E-01 -6.2312E-02  5.3231E-04  2.8266E-01 -1.0017E-01
             2.5134E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1634.30700352477        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1349
 NPARAMETR:  9.9652E-01  1.7424E+00  7.2898E-01  5.2788E-01  1.3264E+00  1.0720E+00  8.4175E-01  1.0000E-02  1.8603E+00  1.0980E+00
             1.0921E+00
 PARAMETER:  9.6567E-02  6.5660E-01 -2.1397E-01 -5.4206E-01  3.8135E-01  1.6959E-01 -7.1373E-02 -4.9412E+00  7.1809E-01  1.9175E-01
             1.8714E-01
 GRADIENT:   3.0503E-02  7.4576E-01  1.1271E-01 -3.5572E-01 -3.6305E-01  5.3433E-03  3.7963E-02  0.0000E+00 -1.3772E-01 -7.7175E-02
            -1.2820E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1349
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.2502E-05 -3.7443E-02 -2.3726E-04  2.8292E-02 -4.4012E-02
 SE:             2.9831E-02  2.1970E-02  9.0929E-05  2.3130E-02  2.2181E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9886E-01  8.8324E-02  9.0720E-03  2.2125E-01  4.7234E-02

 ETASHRINKSD(%)  6.2965E-02  2.6398E+01  9.9695E+01  2.2513E+01  2.5689E+01
 ETASHRINKVR(%)  1.2589E-01  4.5827E+01  9.9999E+01  3.9958E+01  4.4779E+01
 EBVSHRINKSD(%)  4.5993E-01  2.4478E+01  9.9722E+01  2.4519E+01  2.4295E+01
 EBVSHRINKVR(%)  9.1774E-01  4.2965E+01  9.9999E+01  4.3026E+01  4.2688E+01
 RELATIVEINF(%)  9.9016E+01  5.5864E+00  1.8660E-04  6.4148E+00  1.9209E+01
 EPSSHRINKSD(%)  4.2798E+01
 EPSSHRINKVR(%)  6.7279E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1634.3070035247724     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -899.15617696103425     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.70
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1634.307       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.97E-01  1.74E+00  7.31E-01  5.26E-01  1.32E+00  1.07E+00  8.43E-01  1.00E-02  1.86E+00  1.10E+00  1.09E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        9.65E+02
 
 TH 2
+       -5.08E+00  3.17E+02
 
 TH 3
+        5.15E+00  9.76E+01  1.68E+02
 
 TH 4
+       -9.27E+00  2.68E+02 -1.19E+02  6.86E+02
 
 TH 5
+       -3.05E+00 -1.14E+02 -1.29E+02  1.16E+02  3.05E+02
 
 TH 6
+        3.68E-01 -5.90E-01  1.15E+00 -3.00E+00 -7.43E-01  1.70E+02
 
 TH 7
+        7.80E-01  7.31E+00  1.28E+01 -2.22E+01 -1.65E+01 -3.03E-01  1.01E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.35E-01 -1.77E+01 -1.90E+01  4.48E+01  1.24E+00 -3.78E-01  1.32E+01  0.00E+00  2.58E+01
 
 TH10
+        2.65E-01 -9.30E+00 -1.87E+01 -4.32E+00 -4.77E+01  3.85E-02  8.40E+00  0.00E+00  3.29E+00  6.25E+01
 
 TH11
+       -7.68E+00 -2.15E+01 -2.70E+01  3.26E+00 -1.11E+00  2.76E+00  1.10E+01  0.00E+00  2.58E+00  2.01E+01  1.84E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       26.745
Stop Time:
Wed Sep 29 11:12:06 CDT 2021
