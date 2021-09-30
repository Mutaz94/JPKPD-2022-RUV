Wed Sep 29 20:56:42 CDT 2021
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
$DATA ../../../../data/spa1/B/dat21.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2165.51797071118        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6025E+02 -2.2318E+00 -7.3083E+00  2.4588E+01 -3.5405E+01  5.4346E+01  2.9975E+00  1.5860E+01  2.2266E+01  2.0214E+01
            -3.1239E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2176.07377390878        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0471E+00  1.1118E+00  1.1575E+00  9.9914E-01  1.1504E+00  9.9487E-01  1.0193E+00  8.9443E-01  9.0747E-01  9.2601E-01
             1.0303E+00
 PARAMETER:  1.4599E-01  2.0598E-01  2.4622E-01  9.9137E-02  2.4014E-01  9.4856E-02  1.1909E-01 -1.1574E-02  2.9089E-03  2.3132E-02
             1.2989E-01
 GRADIENT:  -3.1779E+00  1.9221E+01  9.1080E-01  2.6635E+01  3.3681E+00 -5.0287E-01 -5.7861E+00  3.3926E+00 -3.7329E+00 -1.6139E+01
            -1.9630E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2177.70173814132        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0494E+00  1.0274E+00  1.2614E+00  1.0377E+00  1.1982E+00  9.9639E-01  1.0615E+00  6.8907E-01  9.1249E-01  1.0495E+00
             1.0277E+00
 PARAMETER:  1.4820E-01  1.2706E-01  3.3219E-01  1.3704E-01  2.8085E-01  9.6385E-02  1.5971E-01 -2.7241E-01  8.4211E-03  1.4829E-01
             1.2736E-01
 GRADIENT:   5.3703E+00 -7.3537E+00 -5.8343E+00  3.1251E+00  2.5822E+01  6.9293E-01 -3.1746E+00 -3.9436E-01  1.0290E+00 -9.8120E+00
            -2.2499E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2178.63135127279        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0474E+00  9.8127E-01  1.1735E+00  1.0617E+00  1.1185E+00  9.9437E-01  1.2028E+00  5.5536E-01  8.5010E-01  1.0458E+00
             1.0516E+00
 PARAMETER:  1.4635E-01  8.1096E-02  2.5997E-01  1.5983E-01  2.1203E-01  9.4353E-02  2.8461E-01 -4.8814E-01 -6.2398E-02  1.4482E-01
             1.5027E-01
 GRADIENT:   7.3336E-01 -1.3438E+00 -1.2685E+00 -6.3086E-01  9.3843E-01  4.5088E-02  3.2645E-02  3.4687E-01  2.6317E-01  3.6507E-01
             3.2431E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2178.78901195013        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0487E+00  1.1567E+00  9.9586E-01  9.4824E-01  1.1156E+00  9.9540E-01  1.0856E+00  2.8298E-01  9.0133E-01  1.0133E+00
             1.0482E+00
 PARAMETER:  1.4751E-01  2.4558E-01  9.5851E-02  4.6848E-02  2.0942E-01  9.5392E-02  1.8212E-01 -1.1624E+00 -3.8825E-03  1.1317E-01
             1.4707E-01
 GRADIENT:  -1.3736E+00  3.2058E+00  1.5773E+00  2.1765E+00 -1.3566E+00 -5.5302E-01  2.9131E-01  6.1230E-02 -7.3685E-01 -8.6807E-01
            -1.4516E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2178.87121097001        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      902             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0529E+00  1.2647E+00  9.1306E-01  8.7588E-01  1.1318E+00  9.9843E-01  1.0118E+00  6.0245E-02  9.5220E-01  1.0097E+00
             1.0479E+00
 PARAMETER:  1.5154E-01  3.3483E-01  9.0442E-03 -3.2521E-02  2.2385E-01  9.8432E-02  1.1172E-01 -2.7093E+00  5.1023E-02  1.0965E-01
             1.4682E-01
 GRADIENT:   6.5356E+02  1.8757E+02  2.8501E+00  5.8391E+01  1.0579E+01  4.6736E+01  5.0999E+00  1.6214E-02  4.2353E+00  1.4027E-01
             8.9479E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2178.87669214474        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  1.0523E+00  1.2672E+00  9.0940E-01  8.7518E-01  1.1319E+00  9.9886E-01  1.0118E+00  3.7887E-02  9.5292E-01  1.0121E+00
             1.0482E+00
 PARAMETER:  1.5097E-01  3.3679E-01  5.0260E-03 -3.3328E-02  2.2389E-01  9.8863E-02  1.1174E-01 -3.1731E+00  5.1776E-02  1.1201E-01
             1.4711E-01
 GRADIENT:   4.1159E+00  1.1187E+00  6.0617E-01  3.5507E-01 -6.6989E-01  3.2389E-01  5.8950E-02  3.2676E-03 -1.0010E-01 -7.6941E-03
            -3.8550E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2178.87782642441        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1258
 NPARAMETR:  1.0522E+00  1.2682E+00  9.0502E-01  8.7320E-01  1.1321E+00  9.9847E-01  1.0101E+00  1.0000E-02  9.5502E-01  1.0113E+00
             1.0482E+00
 PARAMETER:  1.5085E-01  3.3760E-01  2.0368E-04 -3.5594E-02  2.2412E-01  9.8469E-02  1.1000E-01 -4.6078E+00  5.3978E-02  1.1120E-01
             1.4707E-01
 GRADIENT:   3.7751E+00 -1.3637E+00 -7.3336E-01 -3.4275E-01  1.4390E+00  1.5264E-01  1.0631E-02  0.0000E+00  1.6262E-01  9.4846E-02
             1.7978E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2178.87838204482        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  1.0516E+00  1.2686E+00  9.0593E-01  8.7326E-01  1.1318E+00  9.9835E-01  1.0100E+00  1.0000E-02  9.5437E-01  1.0111E+00
             1.0481E+00
 PARAMETER:  1.5033E-01  3.3793E-01  1.2037E-03 -3.5521E-02  2.2379E-01  9.8350E-02  1.0998E-01 -4.6078E+00  5.3295E-02  1.1103E-01
             1.4700E-01
 GRADIENT:   2.5976E+00 -4.1149E-01  2.4146E-02 -4.5676E-01  2.2163E-01  1.0773E-01 -6.1790E-03  0.0000E+00 -6.9763E-03  2.2274E-02
             2.7382E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1315
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.3809E-04 -1.1200E-02 -3.5740E-04  7.0785E-03 -2.5653E-02
 SE:             2.9864E-02  2.1868E-02  1.4398E-04  2.3341E-02  2.2919E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7761E-01  6.0854E-01  1.3052E-02  7.6169E-01  2.6302E-01

 ETASHRINKSD(%)  1.0000E-10  2.6740E+01  9.9518E+01  2.1804E+01  2.3218E+01
 ETASHRINKVR(%)  1.0000E-10  4.6330E+01  9.9998E+01  3.8853E+01  4.1046E+01
 EBVSHRINKSD(%)  3.6725E-01  2.5826E+01  9.9564E+01  2.3302E+01  2.1155E+01
 EBVSHRINKVR(%)  7.3315E-01  4.4982E+01  9.9998E+01  4.1174E+01  3.7835E+01
 RELATIVEINF(%)  9.8737E+01  2.3689E+00  1.9788E-04  2.5574E+00  9.5376E+00
 EPSSHRINKSD(%)  3.2226E+01
 EPSSHRINKVR(%)  5.4067E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2178.8783820448193     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1259.9398488401466     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.21
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2178.878       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.27E+00  9.06E-01  8.73E-01  1.13E+00  9.98E-01  1.01E+00  1.00E-02  9.54E-01  1.01E+00  1.05E+00
 


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
+        1.00E+03
 
 TH 2
+       -6.84E+00  3.87E+02
 
 TH 3
+        7.28E+00  1.13E+02  2.74E+02
 
 TH 4
+       -4.00E+00  3.98E+02 -2.13E+02  9.44E+02
 
 TH 5
+        1.28E+00 -1.80E+02 -3.10E+02  2.31E+02  5.26E+02
 
 TH 6
+        1.73E+00 -1.47E+00  2.15E+00 -1.81E+00 -2.80E-01  1.97E+02
 
 TH 7
+        1.19E+00  2.06E+01  4.31E+00 -1.10E+01 -1.17E+01 -1.09E-01  6.40E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.63E-01 -1.61E+01 -2.57E+01  2.91E+01  1.37E+01 -6.82E-01  3.35E+01  0.00E+00  8.03E+01
 
 TH10
+        1.28E+00 -9.37E+00 -3.19E+01  7.89E-01 -5.24E+01  1.72E-01  1.11E+00  0.00E+00  9.50E+00  8.28E+01
 
 TH11
+       -8.25E+00 -1.80E+01 -3.64E+01 -2.74E+00  3.66E+00  2.26E+00  7.33E+00  0.00E+00  7.27E+00  2.17E+01  3.80E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.803
Stop Time:
Wed Sep 29 20:57:12 CDT 2021
