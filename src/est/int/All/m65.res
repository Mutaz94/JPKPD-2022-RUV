Wed Sep 29 10:41:41 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/28/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/int/All/dat65.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK

ET1 = EXP(ETA(1)*THETA(4))
ET2 = EXP(ETA(2)*THETA(5))
ET3 = EXP(ETA(3)*THETA(6))


CL = 5.0 * THETA(1) * ET1
V = 85  * THETA(2) * ET2
KA = 0.7 * THETA(3) * ET3

SC = V
$ERROR
CVERR 	= 0.05
W  	= THETA(7)*F*CVERR
Y  	= F + W * ERR(1)
$THETA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvKA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvK
(0,1) ; RUV
$OMEGA (0.09 FIX)x3
$SIGMA  1  FIX;        [P]
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
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
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
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
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

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   30874.6086542603        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   6.0000E+02  1.0768E+02 -1.3817E+03 -9.7121E+02 -2.2334E+03 -1.4270E+03 -6.3228E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1090.25991293121        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:       59
 NPARAMETR:  9.8948E-01  2.6458E+00  6.5012E+00  4.1029E+00  5.9278E+00  1.7747E+00  1.0646E+01
 PARAMETER:  8.9421E-02  1.0730E+00  1.9720E+00  1.5117E+00  1.8796E+00  6.7366E-01  2.4652E+00
 GRADIENT:  -1.1448E+01  6.4789E+01  5.8830E+01  1.7422E+02  2.3577E+02  1.8048E+01  2.8546E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1147.71898846196        NO. OF FUNC. EVALS.:  90
 CUMULATIVE NO. OF FUNC. EVALS.:      149             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0066E+00  3.1863E+00  1.9292E+01  1.8111E+00  3.5901E+00  5.6672E+00  9.7448E+00
 PARAMETER:  1.0658E-01  1.2589E+00  3.0597E+00  6.9396E-01  1.3782E+00  1.8347E+00  2.3767E+00
 GRADIENT:  -2.2095E+01  1.4770E+02  2.7712E+00 -3.9230E-02  4.9411E+01 -1.4693E+00  7.7087E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1165.41575653613        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      201
 NPARAMETR:  9.4937E-01  1.7953E+00  1.9901E+01  1.5539E+00  3.4434E+00  5.1118E+00  8.9645E+00
 PARAMETER:  4.8046E-02  6.8516E-01  3.0908E+00  5.4075E-01  1.3364E+00  1.7316E+00  2.2933E+00
 GRADIENT:  -2.1796E+01  3.3838E+01  1.6301E+00 -1.4925E+01  5.8952E+01  1.9058E-01 -9.0022E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1170.18610885872        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      293
 NPARAMETR:  9.9014E-01  1.8297E+00  1.9834E+01  1.5326E+00  3.4546E+00  5.1583E+00  9.5753E+00
 PARAMETER:  9.0091E-02  7.0415E-01  3.0874E+00  5.2698E-01  1.3397E+00  1.7406E+00  2.3592E+00
 GRADIENT:  -1.2583E+01 -4.3097E-01  1.6352E+00 -2.3229E+01 -2.1066E+01  1.7048E-01  2.1754E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1171.02171364702        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:      404
 NPARAMETR:  9.9333E-01  1.8658E+00  1.8372E+01  1.5688E+00  3.4995E+00  5.1793E+00  9.5948E+00
 PARAMETER:  9.3308E-02  7.2368E-01  3.0108E+00  5.5031E-01  1.3526E+00  1.7447E+00  2.3612E+00
 GRADIENT:  -1.1098E+01  2.9020E+00 -4.3528E-01 -1.5913E+01 -1.5645E+01  1.5056E+00  5.3223E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1172.18730015167        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  9.9387E-01  1.8432E+00  1.7929E+01  1.5807E+00  3.5711E+00  5.2275E+00  9.5266E+00
 PARAMETER:  9.3853E-02  7.1149E-01  2.9864E+00  5.5787E-01  1.3729E+00  1.7539E+00  2.3541E+00
 GRADIENT:  -6.7258E+00  1.4852E+00  3.5086E+00 -8.5765E+00 -8.1063E+00 -7.9704E+00  3.5646E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1173.46088889642        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      625
 NPARAMETR:  1.0075E+00  1.7928E+00  1.8046E+01  1.5712E+00  3.6685E+00  5.2971E+00  9.5289E+00
 PARAMETER:  1.0744E-01  6.8379E-01  2.9929E+00  5.5185E-01  1.3998E+00  1.7672E+00  2.3543E+00
 GRADIENT:   8.3495E+00  2.8380E+01  2.2722E+00 -3.1198E+00  8.5584E+01 -3.3634E+00  4.2486E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1173.51008462803        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:      758
 NPARAMETR:  1.0047E+00  1.8299E+00  1.8119E+01  1.5707E+00  3.6692E+00  5.3511E+00  9.5691E+00
 PARAMETER:  1.0465E-01  7.0428E-01  2.9970E+00  5.5151E-01  1.4000E+00  1.7773E+00  2.3585E+00
 GRADIENT:  -2.0308E+00  5.9018E-03  1.3209E+00 -1.3190E+01  2.7145E+00 -5.2810E+00  1.1085E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1173.60965223868        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      875
 NPARAMETR:  9.9870E-01  1.8424E+00  1.8176E+01  1.5850E+00  3.6698E+00  5.3704E+00  9.5456E+00
 PARAMETER:  9.8703E-02  7.1107E-01  3.0001E+00  5.6060E-01  1.4001E+00  1.7809E+00  2.3561E+00
 GRADIENT:   5.7903E-03  1.1146E+00  1.5659E+00 -2.0587E+00  2.2916E+00 -5.4028E+00  1.8360E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1173.66058166173        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1010
 NPARAMETR:  9.9274E-01  1.8366E+00  1.8205E+01  1.5933E+00  3.6749E+00  5.4289E+00  9.5555E+00
 PARAMETER:  9.2713E-02  7.0790E-01  3.0017E+00  5.6579E-01  1.4015E+00  1.7917E+00  2.3571E+00
 GRADIENT:  -2.2920E+00  5.2475E-01  7.1195E-01  1.7937E+00  2.6142E+00 -4.8505E+00  3.8253E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1173.69737945657        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:     1136
 NPARAMETR:  9.9815E-01  1.8395E+00  1.8330E+01  1.5874E+00  3.6593E+00  5.4844E+00  9.5827E+00
 PARAMETER:  9.9149E-02  7.1097E-01  3.0042E+00  5.6247E-01  1.3941E+00  1.8052E+00  2.3571E+00
 GRADIENT:   5.5684E+00  2.9459E+01 -9.2801E+01  6.6832E+01 -1.7767E+01  3.5239E+01 -4.6110E+01
 NUMSIGDIG:         1.6         2.3         2.5         2.8         2.3         2.4         2.5

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1136
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5866E-02 -8.4369E-03 -2.5717E-02
 SE:             2.8625E-02  3.0031E-02  7.1250E-03
 N:                     100         100         100

 P VAL.:         3.6619E-01  7.7875E-01  3.0701E-04

 ETASHRINKSD(%)  4.1035E+00  1.0000E-10  7.6130E+01
 ETASHRINKVR(%)  8.0387E+00  1.0000E-10  9.4302E+01
 EBVSHRINKSD(%)  1.3244E+01  1.1115E+00  8.4577E+01
 EBVSHRINKVR(%)  2.4734E+01  2.2107E+00  9.7621E+01
 RELATIVEINF(%)  7.4906E+01  9.7779E+01  2.3674E+00
 EPSSHRINKSD(%)  7.9368E+00
 EPSSHRINKVR(%)  1.5244E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1173.6973794565713     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       480.39198031183946     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:    30.32
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1173.697       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         9.99E-01  1.84E+00  1.83E+01  1.59E+00  3.65E+00  5.50E+00  9.55E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 


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
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        1.04E+05
 
 TH 2
+        7.74E+03  6.17E+02
 
 TH 3
+        1.13E+00 -3.31E-01  9.26E+00
 
 TH 4
+        1.04E+04  9.97E+00 -5.56E-01  3.54E+03
 
 TH 5
+       -4.29E+01  5.14E-01  3.53E-01 -2.03E+02  5.49E+01
 
 TH 6
+       -5.42E-01  3.94E-01 -7.91E+00  1.25E-01 -4.68E-01  6.94E+01
 
 TH 7
+        1.43E+01  4.22E+01 -3.48E-01  8.97E+01  1.14E+01 -1.38E+01  3.23E+01
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       36.774
Stop Time:
Wed Sep 29 10:42:19 CDT 2021
