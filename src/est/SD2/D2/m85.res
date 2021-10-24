Sat Oct 23 20:50:37 CDT 2021
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
$DATA ../../../../data/SD2/D2/dat85.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2894.04968955462        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5294E+02 -2.2868E+02 -5.6972E+01 -1.7001E+02  1.1774E+02  1.0318E+01 -1.2325E+02  1.8339E+01 -1.5155E+02  1.2335E+01
            -4.7171E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2918.43923648176        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:       96
 NPARAMETR:  9.2086E-01  1.1069E+00  1.1772E+00  9.3855E-01  1.0692E+00  9.4175E-01  1.6899E+00  8.1915E-01  9.0446E-01  7.7456E-01
             9.7129E-01
 PARAMETER:  1.7556E-02  2.0157E-01  2.6314E-01  3.6583E-02  1.6693E-01  3.9986E-02  6.2464E-01 -9.9491E-02 -4.1517E-04 -1.5546E-01
             7.0874E-02
 GRADIENT:   5.4748E+01  2.8310E+01  7.0375E+01 -1.4862E+02 -7.1822E+01 -4.7540E+01  9.0661E+01 -8.7313E-01 -9.5616E+01 -1.6294E+01
            -1.1006E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2954.80660304832        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  1.1112E+00  1.4738E+00  9.6244E-01  1.0162E+00  1.1187E+00  8.6517E-01  1.8554E+00  7.3519E-01  8.9183E-01  8.6115E-01
             1.0161E+00
 PARAMETER:  2.0546E-01  4.8784E-01  6.1719E-02  1.1602E-01  2.1215E-01 -4.4834E-02  7.1811E-01 -2.0763E-01 -1.4483E-02 -4.9482E-02
             1.1597E-01
 GRADIENT:   9.5072E+01  5.8876E+01 -6.4527E+01  1.1235E+02  6.2643E+01 -1.2729E+02  1.6477E+01  5.6115E+00 -6.7785E+01 -1.9032E+01
            -4.1898E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2994.23351513229        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      415
 NPARAMETR:  1.0663E+00  1.3382E+00  1.0486E+00  9.7043E-01  1.0959E+00  1.0793E+00  1.5062E+00  3.2536E-01  1.2373E+00  8.3007E-01
             1.0410E+00
 PARAMETER:  1.6421E-01  3.9133E-01  1.4747E-01  6.9989E-02  1.9160E-01  1.7630E-01  5.0957E-01 -1.0228E+00  3.1291E-01 -8.6241E-02
             1.4014E-01
 GRADIENT:  -1.6279E+01 -5.8329E+00  1.7370E+01 -1.0010E+01 -2.4235E+01 -9.5552E+00 -1.7729E+01 -5.9563E-01 -1.9930E+01 -5.1827E+00
             4.0617E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2996.55632547379        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.0766E+00  1.4487E+00  1.0261E+00  9.1701E-01  1.1580E+00  1.1029E+00  1.5232E+00  3.4285E-01  1.3813E+00  8.7876E-01
             1.0408E+00
 PARAMETER:  1.7377E-01  4.7067E-01  1.2580E-01  1.3366E-02  2.4666E-01  1.9795E-01  5.2081E-01 -9.7047E-01  4.2305E-01 -2.9242E-02
             1.4003E-01
 GRADIENT:   2.3487E+00  9.0822E-01 -2.0105E+00  4.0333E+00  9.0198E-01 -5.5983E-01 -5.4832E-01  1.9537E-01  4.3405E-03 -3.4255E-02
             1.7524E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2996.57969129354        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      767             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0778E+00  1.4580E+00  1.0243E+00  9.0661E-01  1.1648E+00  1.1061E+00  1.5224E+00  3.1689E-01  1.3906E+00  8.8245E-01
             1.0400E+00
 PARAMETER:  1.7489E-01  4.7708E-01  1.2398E-01  1.9545E-03  2.5251E-01  2.0080E-01  5.2028E-01 -1.0492E+00  4.2971E-01 -2.5057E-02
             1.3920E-01
 GRADIENT:   6.8269E+02  3.5936E+02  2.9885E+00  5.7487E+01  4.8907E+01  1.0509E+02  7.5453E+01  2.7646E-01  3.5243E+01  1.3243E+00
             1.9445E+00

0ITERATION NO.:   27    OBJECTIVE VALUE:  -2996.57969129354        NO. OF FUNC. EVALS.:  62
 CUMULATIVE NO. OF FUNC. EVALS.:      829
 NPARAMETR:  1.0778E+00  1.4580E+00  1.0255E+00  9.0619E-01  1.1638E+00  1.1061E+00  1.5227E+00  3.1358E-01  1.3899E+00  8.8270E-01
             1.0400E+00
 PARAMETER:  1.7489E-01  4.7708E-01  1.2398E-01  1.9545E-03  2.5251E-01  2.0080E-01  5.2028E-01 -1.0492E+00  4.2971E-01 -2.5057E-02
             1.3920E-01
 GRADIENT:  -5.2085E-03  1.1898E-02 -3.4827E-01  3.9033E-01  1.3170E+00 -1.0825E-03 -5.0534E-02  1.2209E-02  1.0244E-01 -4.0497E-02
            -3.4433E-03
 NUMSIGDIG:         4.8         4.5         1.9         2.3         2.5         4.9         3.4         1.9         2.9         2.5
                    4.8

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      829
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5431E-04 -1.2078E-02 -8.9607E-03  1.2623E-02 -2.0627E-02
 SE:             2.9914E-02  2.5898E-02  5.0444E-03  2.4839E-02  2.1989E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9322E-01  6.4094E-01  7.5671E-02  6.1133E-01  3.4820E-01

 ETASHRINKSD(%)  1.0000E-10  1.3237E+01  8.3101E+01  1.6785E+01  2.6335E+01
 ETASHRINKVR(%)  1.0000E-10  2.4722E+01  9.7144E+01  3.0753E+01  4.5735E+01
 EBVSHRINKSD(%)  2.4055E-01  1.3290E+01  8.4652E+01  1.7080E+01  2.6611E+01
 EBVSHRINKVR(%)  4.8053E-01  2.4814E+01  9.7644E+01  3.1243E+01  4.6141E+01
 RELATIVEINF(%)  9.9517E+01  2.6668E+01  1.4912E+00  2.3679E+01  1.6883E+01
 EPSSHRINKSD(%)  2.4237E+01
 EPSSHRINKVR(%)  4.2600E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2996.5796912935371     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1710.0657448069953     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2996.580       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.08E+00  1.46E+00  1.02E+00  9.07E-01  1.16E+00  1.11E+00  1.52E+00  3.17E-01  1.39E+00  8.82E-01  1.04E+00
 


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
 #CPUT: Total CPU Time in Seconds,       73.789
Stop Time:
Sat Oct 23 20:50:49 CDT 2021
