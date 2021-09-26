Sat Sep 25 10:38:18 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat63.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1693.37313869068        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.8096E+01 -1.1233E+02 -4.5023E+01 -1.2929E+02  3.3278E+01 -2.2934E+01 -1.1780E+01  1.3784E+01 -1.4972E+01  2.3296E+01
             1.4922E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1708.23553906349        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      133
 NPARAMETR:  1.0508E+00  1.0906E+00  1.1595E+00  1.0339E+00  1.0667E+00  1.0869E+00  1.0444E+00  9.2690E-01  1.0385E+00  8.7339E-01
             9.5305E-01
 PARAMETER:  1.4954E-01  1.8668E-01  2.4801E-01  1.3334E-01  1.6454E-01  1.8331E-01  1.4342E-01  2.4092E-02  1.3781E-01 -3.5374E-02
             5.1908E-02
 GRADIENT:   7.8087E+00 -1.1452E+01  6.3780E+00 -2.0347E+01  9.9983E+00  9.7805E+00 -4.2765E+00 -1.3835E+00 -3.2319E+00 -8.6537E+00
            -1.3207E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1708.83735700354        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      310
 NPARAMETR:  1.0485E+00  1.1062E+00  1.0634E+00  1.0267E+00  1.0366E+00  1.0713E+00  1.1943E+00  7.5463E-01  9.9346E-01  8.4720E-01
             9.7094E-01
 PARAMETER:  1.4735E-01  2.0091E-01  1.6149E-01  1.2634E-01  1.3591E-01  1.6884E-01  2.7759E-01 -1.8153E-01  9.3439E-02 -6.5822E-02
             7.0512E-02
 GRADIENT:   2.5835E+00 -9.9306E-01  6.7628E+00 -1.2610E+01  1.3068E+01  3.9710E+00  3.9812E+00 -1.4155E+00 -4.4238E+00 -6.3858E+00
            -5.3444E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1710.23053626182        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  1.0477E+00  9.3002E-01  8.9550E-01  1.1325E+00  8.6873E-01  1.0582E+00  1.3571E+00  4.7017E-01  9.1781E-01  7.5323E-01
             9.8038E-01
 PARAMETER:  1.4663E-01  2.7448E-02 -1.0377E-02  2.2439E-01 -4.0724E-02  1.5654E-01  4.0536E-01 -6.5465E-01  1.4230E-02 -1.8338E-01
             8.0183E-02
 GRADIENT:   1.1875E-01  4.4113E+00 -1.3236E+00  8.6110E+00  3.7767E-01 -9.9230E-01 -3.8226E-01  5.4189E-01  3.3417E-01  5.5517E-01
             9.0081E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1710.43745081471        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  1.0475E+00  8.2174E-01  8.2405E-01  1.1783E+00  7.8955E-01  1.0605E+00  1.5292E+00  3.0074E-01  8.6811E-01  6.9187E-01
             9.8079E-01
 PARAMETER:  1.4642E-01 -9.6328E-02 -9.3529E-02  2.6404E-01 -1.3630E-01  1.5877E-01  5.2474E-01 -1.1015E+00 -4.1438E-02 -2.6836E-01
             8.0602E-02
 GRADIENT:  -2.0634E-01 -9.7410E-01 -8.8773E-01 -1.0294E+00  7.3564E-01 -1.6239E-01  3.3289E-01  1.3140E-01  2.5985E-01  4.2179E-01
             4.3230E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1710.47694787498        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      842
 NPARAMETR:  1.0480E+00  8.9273E-01  7.8321E-01  1.1332E+00  7.9416E-01  1.0616E+00  1.4395E+00  1.6925E-01  8.8676E-01  6.8186E-01
             9.8090E-01
 PARAMETER:  1.4685E-01 -1.3470E-02 -1.4436E-01  2.2502E-01 -1.3047E-01  1.5976E-01  4.6430E-01 -1.6764E+00 -2.0178E-02 -2.8292E-01
             8.0713E-02
 GRADIENT:  -1.0019E+00 -2.8610E-02 -1.5028E-01 -2.8408E-01  7.6170E-01 -5.0919E-02  3.9330E-01  2.0632E-02  1.7648E-01 -2.0895E-01
             2.3874E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1710.48659833825        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1017
 NPARAMETR:  1.0486E+00  8.9421E-01  7.6815E-01  1.1304E+00  7.8615E-01  1.0617E+00  1.4362E+00  7.6884E-02  8.8466E-01  6.7664E-01
             9.8075E-01
 PARAMETER:  1.4741E-01 -1.1818E-02 -1.6377E-01  2.2256E-01 -1.4061E-01  1.5987E-01  4.6203E-01 -2.4655E+00 -2.2551E-02 -2.9061E-01
             8.0561E-02
 GRADIENT:  -2.3591E-01 -1.5387E-01 -2.5072E-01  9.8419E-03  3.6374E-01 -4.6128E-02  8.4696E-02  3.6386E-03  8.5582E-02 -9.0681E-03
             1.3213E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1710.48833408271        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1192
 NPARAMETR:  1.0487E+00  9.0025E-01  7.6485E-01  1.1266E+00  7.8655E-01  1.0619E+00  1.4282E+00  2.0009E-02  8.8608E-01  6.7659E-01
             9.8050E-01
 PARAMETER:  1.4757E-01 -5.0868E-03 -1.6808E-01  2.1923E-01 -1.4009E-01  1.6002E-01  4.5644E-01 -3.8116E+00 -2.0947E-02 -2.9069E-01
             8.0308E-02
 GRADIENT:  -4.0740E-02  2.3083E-02  2.4299E-03  3.0123E-02 -7.1934E-03 -5.5088E-03  1.0623E-02  2.1666E-04  5.8673E-03 -3.3018E-03
             1.4726E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1710.48841861509        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1319
 NPARAMETR:  1.0487E+00  9.0050E-01  7.6468E-01  1.1264E+00  7.8657E-01  1.0619E+00  1.4278E+00  1.0000E-02  8.8614E-01  6.7662E-01
             9.8046E-01
 PARAMETER:  1.4760E-01 -4.8039E-03 -1.6830E-01  2.1906E-01 -1.4007E-01  1.6004E-01  4.5614E-01 -4.7124E+00 -2.0879E-02 -2.9064E-01
             8.0271E-02
 GRADIENT:  -7.4365E-04 -3.8337E-03 -2.7165E-03 -1.5911E-03  5.6629E-03 -1.6567E-04 -1.1897E-04  0.0000E+00  3.7902E-04  2.4851E-04
            -1.7122E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1319
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0465E-04  7.9407E-03 -5.2227E-04 -8.6680E-03 -4.7461E-03
 SE:             2.9856E-02  2.2493E-02  2.2133E-04  2.5026E-02  2.0999E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9720E-01  7.2407E-01  1.8288E-02  7.2907E-01  8.2119E-01

 ETASHRINKSD(%)  1.0000E-10  2.4644E+01  9.9259E+01  1.6161E+01  2.9652E+01
 ETASHRINKVR(%)  1.0000E-10  4.3215E+01  9.9995E+01  2.9710E+01  5.0511E+01
 EBVSHRINKSD(%)  3.7168E-01  2.4271E+01  9.9318E+01  1.6106E+01  2.9359E+01
 EBVSHRINKVR(%)  7.4199E-01  4.2652E+01  9.9995E+01  2.9618E+01  5.0099E+01
 RELATIVEINF(%)  9.9018E+01  5.4906E+00  5.1067E-04  8.7086E+00  3.8256E+00
 EPSSHRINKSD(%)  4.3525E+01
 EPSSHRINKVR(%)  6.8105E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1710.4884186150857     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -975.33759205134754     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.99
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1710.488       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  9.01E-01  7.65E-01  1.13E+00  7.87E-01  1.06E+00  1.43E+00  1.00E-02  8.86E-01  6.77E-01  9.80E-01
 


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
+        8.91E+02
 
 TH 2
+       -6.17E+00  4.02E+02
 
 TH 3
+        1.60E+01  2.96E+02  9.51E+02
 
 TH 4
+       -7.50E+00  2.47E+02 -3.64E+02  7.98E+02
 
 TH 5
+       -6.08E+00 -5.11E+02 -1.21E+03  4.45E+02  1.94E+03
 
 TH 6
+        3.36E-01 -2.00E+00  1.47E+00 -2.36E+00 -1.73E+00  1.74E+02
 
 TH 7
+        8.59E-01  3.27E+01 -1.30E+01 -1.15E+01  2.96E+00 -3.02E-01  3.73E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.54E-01 -2.11E+01 -3.95E+01  2.89E+01 -2.90E+00 -1.47E+00  1.24E+01  0.00E+00  1.38E+02
 
 TH10
+       -1.97E+00 -1.01E+01 -8.73E+01 -3.36E+01 -4.15E+01  1.23E+00  1.57E+01  0.00E+00  1.39E+01  1.19E+02
 
 TH11
+       -4.65E+00 -1.26E+01 -4.49E+01 -7.94E+00  1.14E+01  2.96E+00  4.38E+00  0.00E+00  1.13E+01  2.97E+01  2.28E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.710
Stop Time:
Sat Sep 25 10:38:41 CDT 2021
