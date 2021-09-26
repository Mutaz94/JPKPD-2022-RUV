Sat Sep 25 10:51:22 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat96.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1636.80916147691        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2344E+02 -8.2644E+00 -3.5177E+01  3.6110E+01  3.0898E+01 -9.8853E+00  9.1322E+00  1.0068E+01  2.4768E+01  1.4215E+01
            -4.6028E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1643.68902714973        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.8692E-01  1.0459E+00  1.1107E+00  9.8020E-01  1.0430E+00  1.0239E+00  9.2923E-01  9.2364E-01  8.6189E-01  8.9993E-01
             1.1352E+00
 PARAMETER:  8.6834E-02  1.4484E-01  2.0500E-01  8.0006E-02  1.4206E-01  1.2364E-01  2.6599E-02  2.0569E-02 -4.8628E-02 -5.4354E-03
             2.2677E-01
 GRADIENT:   8.1300E+01  2.3321E+01 -1.2622E+00  3.2539E+01  1.6382E+01  2.3188E+00 -2.5009E+00  8.7855E-01 -8.3560E+00 -8.7106E+00
             2.6815E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1644.51790453915        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.7959E-01  9.9402E-01  1.1187E+00  1.0065E+00  1.0363E+00  1.0324E+00  9.1788E-01  5.4003E-01  9.2957E-01  9.8726E-01
             1.1349E+00
 PARAMETER:  7.9380E-02  9.4001E-02  2.1220E-01  1.0647E-01  1.3570E-01  1.3190E-01  1.4308E-02 -5.1613E-01  2.6967E-02  8.7179E-02
             2.2658E-01
 GRADIENT:   6.6184E+01  1.0312E+01  3.1052E+00  2.3955E+01  1.5137E+01  6.6103E+00  1.2115E+00 -1.2852E+00  6.1842E+00 -2.7337E+00
             2.3329E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1645.68662845709        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.5334E-01  1.0014E+00  8.8837E-01  9.7834E-01  9.2629E-01  1.0137E+00  1.0221E+00  3.7605E-01  8.7771E-01  8.8531E-01
             1.1128E+00
 PARAMETER:  5.2211E-02  1.0140E-01 -1.8366E-02  7.8106E-02  2.3427E-02  1.1359E-01  1.2184E-01 -8.7802E-01 -3.0439E-02 -2.1823E-02
             2.0685E-01
 GRADIENT:   7.2299E+00  1.3167E-01 -3.6681E+00  3.4477E+00  6.0760E+00 -1.2524E+00  2.5463E+00  4.0509E-01  1.8865E+00  1.1498E+00
             2.2970E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1645.99428527712        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      405
 NPARAMETR:  9.6365E-01  9.1790E-01  8.9454E-01  1.0328E+00  8.8804E-01  1.0260E+00  1.0757E+00  2.9352E-01  8.3963E-01  8.6755E-01
             1.1149E+00
 PARAMETER:  6.2977E-02  1.4337E-02 -1.1442E-02  1.3231E-01 -1.8737E-02  1.2566E-01  1.7297E-01 -1.1258E+00 -7.4794E-02 -4.2084E-02
             2.0876E-01
 GRADIENT:  -2.8904E+00  2.1163E+00 -1.3438E-01  3.3006E+00 -3.7681E-01 -4.6360E-01  2.4300E-01  1.1622E-01 -8.0820E-02 -1.2194E-01
             2.3265E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1646.03214824691        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  9.6454E-01  9.1638E-01  8.6837E-01  1.0297E+00  8.7433E-01  1.0265E+00  1.0809E+00  1.4225E-01  8.3998E-01  8.6169E-01
             1.1154E+00
 PARAMETER:  6.3893E-02  1.2671E-02 -4.1136E-02  1.2931E-01 -3.4297E-02  1.2611E-01  1.7776E-01 -1.8502E+00 -7.4379E-02 -4.8854E-02
             2.0923E-01
 GRADIENT:  -1.4231E+00  3.7767E-02  3.3439E-01 -3.9614E-01 -1.2507E+00 -3.6465E-01  2.4377E-01  2.1937E-02  4.1931E-01  2.2633E-01
             3.6010E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1646.04831711151        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      757
 NPARAMETR:  9.6552E-01  9.7515E-01  8.4850E-01  9.9314E-01  8.9104E-01  1.0278E+00  1.0291E+00  5.7126E-02  8.6163E-01  8.6528E-01
             1.1149E+00
 PARAMETER:  6.4907E-02  7.4833E-02 -6.4289E-02  9.3116E-02 -1.5367E-02  1.2741E-01  1.2873E-01 -2.7625E+00 -4.8933E-02 -4.4702E-02
             2.0874E-01
 GRADIENT:  -1.1652E-01 -1.1622E-01 -4.3115E-01  2.7379E-01  4.7444E-01  1.3966E-02 -4.7461E-02  4.2175E-03  2.6551E-02  6.6243E-02
             5.2321E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1646.05054834369        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      932
 NPARAMETR:  9.6553E-01  9.6779E-01  8.5030E-01  9.9763E-01  8.8853E-01  1.0277E+00  1.0358E+00  1.0000E-02  8.5844E-01  8.6468E-01
             1.1149E+00
 PARAMETER:  6.4921E-02  6.7258E-02 -6.2161E-02  9.7625E-02 -1.8186E-02  1.2737E-01  1.3516E-01 -4.5908E+00 -5.2639E-02 -4.5400E-02
             2.0876E-01
 GRADIENT:  -4.5266E-04 -2.8175E-02 -7.7930E-02  5.2112E-02  1.0614E-01  8.4982E-03 -1.3305E-02  0.0000E+00 -6.1413E-03  4.9274E-03
             1.9765E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1646.05055837017        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  9.6552E-01  9.6657E-01  8.5079E-01  9.9838E-01  8.8819E-01  1.0277E+00  1.0369E+00  1.0000E-02  8.5797E-01  8.6463E-01
             1.1149E+00
 PARAMETER:  6.4913E-02  6.5996E-02 -6.1595E-02  9.8382E-02 -1.8572E-02  1.2734E-01  1.3624E-01 -4.5609E+00 -5.3185E-02 -4.5457E-02
             2.0877E-01
 GRADIENT:  -4.5816E-06  2.0974E-04  4.5501E-04 -1.7120E-04 -3.0680E-04 -3.9944E-05  1.5469E-05  0.0000E+00 -9.8294E-06 -1.9920E-04
            -7.4496E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1024
 NO. OF SIG. DIGITS IN FINAL EST.:  4.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.1169E-04 -5.5062E-03 -3.1794E-04 -2.2422E-03 -1.8400E-02
 SE:             2.9803E-02  1.9691E-02  1.6600E-04  2.4173E-02  2.3201E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9433E-01  7.7976E-01  5.5444E-02  9.2610E-01  4.2775E-01

 ETASHRINKSD(%)  1.5471E-01  3.4034E+01  9.9444E+01  1.9019E+01  2.2273E+01
 ETASHRINKVR(%)  3.0919E-01  5.6484E+01  9.9997E+01  3.4421E+01  3.9585E+01
 EBVSHRINKSD(%)  4.9012E-01  3.3969E+01  9.9474E+01  1.9236E+01  2.1049E+01
 EBVSHRINKVR(%)  9.7784E-01  5.6399E+01  9.9997E+01  3.4771E+01  3.7668E+01
 RELATIVEINF(%)  9.8542E+01  1.4315E+00  2.5634E-04  2.6351E+00  4.9967E+00
 EPSSHRINKSD(%)  4.1909E+01
 EPSSHRINKVR(%)  6.6254E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1646.0505583701661     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -910.89973180642789     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.60
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1646.051       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.66E-01  9.67E-01  8.51E-01  9.98E-01  8.88E-01  1.03E+00  1.04E+00  1.00E-02  8.58E-01  8.65E-01  1.11E+00
 


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
+        1.12E+03
 
 TH 2
+       -8.52E+00  5.00E+02
 
 TH 3
+        1.36E+01  2.06E+02  4.47E+02
 
 TH 4
+       -1.43E+01  4.67E+02 -2.08E+02  9.91E+02
 
 TH 5
+       -5.59E+00 -3.97E+02 -6.39E+02  2.41E+02  1.17E+03
 
 TH 6
+       -1.61E+00 -2.68E+00  2.23E+00 -2.80E+00 -2.41E+00  1.84E+02
 
 TH 7
+        1.04E-01  2.45E+01  6.99E+00 -8.05E+00 -1.35E+01 -1.02E+00  4.06E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.92E+00 -2.11E+01 -2.67E+01  2.64E+01  6.76E+00  4.53E-01  2.54E+01  0.00E+00  1.26E+02
 
 TH10
+       -1.11E-01 -5.03E+00 -4.04E+01 -1.76E+01 -5.97E+01  1.05E+00  1.54E+01  0.00E+00  5.78E+00  1.06E+02
 
 TH11
+       -8.10E+00 -1.73E+01 -3.58E+01 -5.56E+00  9.13E+00  2.66E+00  5.95E+00  0.00E+00  1.13E+01  2.44E+01  1.80E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.839
Stop Time:
Sat Sep 25 10:51:41 CDT 2021
