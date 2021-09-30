Wed Sep 29 16:13:19 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat90.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m90.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1618.25005896283        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6460E+02 -3.9084E+00 -2.9054E+01  4.7160E+01  6.5709E+01  2.9603E+01  1.2952E-01  4.9400E+00 -1.2618E+00 -2.1540E+00
            -4.4902E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1622.79026777349        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8048E-01  1.0180E+00  1.0237E+00  1.0130E+00  9.7830E-01  1.0141E+00  1.0046E+00  9.8669E-01  1.0260E+00  9.9136E-01
             1.1048E+00
 PARAMETER:  8.0291E-02  1.1788E-01  1.2344E-01  1.1292E-01  7.8065E-02  1.1399E-01  1.0459E-01  8.6597E-02  1.2567E-01  9.1325E-02
             1.9966E-01
 GRADIENT:   2.1810E-01 -2.3147E+00 -8.3492E+00  5.3627E+00  9.6394E+00 -7.4669E-02  1.1303E-01  3.3351E+00 -2.1054E+00  1.9431E+00
             3.6573E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1622.97802807443        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.8098E-01  1.0294E+00  1.0085E+00  1.0034E+00  9.7645E-01  1.0112E+00  9.7992E-01  8.6566E-01  1.0368E+00  9.8287E-01
             1.1048E+00
 PARAMETER:  8.0794E-02  1.2894E-01  1.0844E-01  1.0341E-01  7.6167E-02  1.1113E-01  7.9716E-02 -4.4258E-02  1.3611E-01  8.2726E-02
             1.9967E-01
 GRADIENT:   1.0683E+00 -1.6504E+00 -1.9761E+00  2.8397E+00  9.3415E+00 -1.1991E+00 -1.3717E+00 -4.7341E-01 -2.4914E+00 -1.4906E+00
            -8.7099E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1623.23691577600        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.8170E-01  1.1952E+00  8.4474E-01  8.9242E-01  9.6663E-01  1.0171E+00  9.3289E-01  7.3358E-01  1.1253E+00  9.5818E-01
             1.1051E+00
 PARAMETER:  8.1528E-02  2.7833E-01 -6.8723E-02 -1.3823E-02  6.6063E-02  1.1692E-01  3.0533E-02 -2.0981E-01  2.1804E-01  5.7285E-02
             1.9990E-01
 GRADIENT:  -5.6509E-01  4.4133E+00  2.1527E+00  3.0279E+00 -5.8090E+00  4.1755E-01  4.8475E-01  1.8121E-02  5.5219E-01  1.1710E+00
             2.4894E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1623.48682472446        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      706
 NPARAMETR:  9.8365E-01  1.5015E+00  6.2495E-01  6.8988E-01  1.0172E+00  1.0168E+00  8.2619E-01  5.1319E-01  1.3248E+00  9.3729E-01
             1.1060E+00
 PARAMETER:  8.3514E-02  5.0647E-01 -3.7009E-01 -2.7124E-01  1.1709E-01  1.1670E-01 -9.0926E-02 -5.6711E-01  3.8125E-01  3.5242E-02
             2.0073E-01
 GRADIENT:   8.8418E-01  9.9802E+00  2.1146E+00  4.8715E+00 -5.2003E+00 -2.2403E-01  1.9425E-01  2.0919E-01  3.2825E-01 -1.1158E+00
            -2.3983E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1623.60021049993        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      887
 NPARAMETR:  9.8384E-01  1.5435E+00  5.7966E-01  6.5201E-01  1.0263E+00  1.0178E+00  8.1001E-01  3.1000E-01  1.3592E+00  9.3947E-01
             1.1055E+00
 PARAMETER:  8.3708E-02  5.3408E-01 -4.4531E-01 -3.2769E-01  1.2596E-01  1.1764E-01 -1.1071E-01 -1.0712E+00  4.0693E-01  3.7562E-02
             2.0028E-01
 GRADIENT:   1.1583E+00 -5.5851E+00  8.9350E-01 -1.8932E+00  2.0903E+00  1.2409E-01 -8.9812E-01 -2.4642E-03 -7.3610E-01 -1.1075E+00
            -2.2259E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1623.62993375439        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  9.8346E-01  1.5760E+00  5.5316E-01  6.3231E-01  1.0258E+00  1.0180E+00  8.0313E-01  1.8830E-01  1.3777E+00  9.3934E-01
             1.1029E+00
 PARAMETER:  8.3322E-02  5.5488E-01 -4.9212E-01 -3.5837E-01  1.2543E-01  1.1782E-01 -1.1923E-01 -1.5697E+00  4.2040E-01  3.7418E-02
             1.9790E-01
 GRADIENT:   1.8878E-01  1.8069E+00  2.8943E+00 -1.3489E+00 -4.2462E+00  1.5119E-01 -8.3715E-01  4.5096E-03 -1.4953E+00 -4.6674E-01
            -1.1663E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1623.63703439331        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.8359E-01  1.6194E+00  5.1853E-01  6.0277E-01  1.0308E+00  1.0181E+00  7.9121E-01  8.2813E-02  1.4121E+00  9.3295E-01
             1.1021E+00
 PARAMETER:  8.3450E-02  5.8207E-01 -5.5675E-01 -4.0622E-01  1.3031E-01  1.1797E-01 -1.3419E-01 -2.3912E+00  4.4506E-01  3.0597E-02
             1.9724E-01
 GRADIENT:   3.4428E-01  2.4440E+00  2.2373E+00 -1.7057E+00 -5.7660E+00  1.7853E-01 -7.7959E-01  1.0484E-02 -1.6764E+00 -5.5460E-01
            -1.3370E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1623.66514313212        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1421
 NPARAMETR:  9.8323E-01  1.6235E+00  5.1820E-01  6.0059E-01  1.0375E+00  1.0175E+00  7.9255E-01  2.0761E-02  1.4269E+00  9.3980E-01
             1.1054E+00
 PARAMETER:  8.3087E-02  5.8459E-01 -5.5740E-01 -4.0984E-01  1.3679E-01  1.1731E-01 -1.3250E-01 -3.7747E+00  4.5547E-01  3.7910E-02
             2.0018E-01
 GRADIENT:  -4.8708E-01 -1.3125E+00 -8.4265E-02  5.7059E-01 -6.6895E-01 -7.8261E-02  4.5482E-02  7.8059E-04 -1.4279E-01 -7.9136E-02
             1.6041E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1623.66759009174        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:     1525
 NPARAMETR:  9.8345E-01  1.6237E+00  5.2077E-01  6.0123E-01  1.0371E+00  1.0179E+00  7.9202E-01  1.0000E-02  1.4287E+00  9.4020E-01
             1.1044E+00
 PARAMETER:  8.4307E-02  5.8332E-01 -5.5710E-01 -4.1105E-01  1.3725E-01  1.1778E-01 -1.3269E-01 -4.6259E+00  4.5669E-01  3.8512E-02
             1.9980E-01
 GRADIENT:   8.2382E-01 -1.1757E+00 -4.4885E-01 -5.9814E-01  4.3525E-01  8.0307E-03  2.9284E-02  0.0000E+00 -6.8671E-03  1.0035E-02
             8.5410E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1525
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.6896E-04 -2.9880E-02 -2.8025E-04  2.3449E-02 -3.5003E-02
 SE:             2.9822E-02  2.2895E-02  1.1171E-04  2.3685E-02  2.2078E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9548E-01  1.9186E-01  1.2118E-02  3.2215E-01  1.1286E-01

 ETASHRINKSD(%)  9.4208E-02  2.3298E+01  9.9626E+01  2.0653E+01  2.6037E+01
 ETASHRINKVR(%)  1.8833E-01  4.1168E+01  9.9999E+01  3.7040E+01  4.5295E+01
 EBVSHRINKSD(%)  5.0487E-01  2.2820E+01  9.9677E+01  2.1710E+01  2.4985E+01
 EBVSHRINKVR(%)  1.0072E+00  4.0433E+01  9.9999E+01  3.8706E+01  4.3727E+01
 RELATIVEINF(%)  9.8961E+01  4.1348E+00  1.2907E-04  4.7003E+00  1.0043E+01
 EPSSHRINKSD(%)  4.4000E+01
 EPSSHRINKVR(%)  6.8640E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1623.6675900917435     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -888.51676352800530     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.70
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1623.668       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.62E+00  5.18E-01  6.00E-01  1.04E+00  1.02E+00  7.92E-01  1.00E-02  1.43E+00  9.40E-01  1.10E+00
 


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
+        1.10E+03
 
 TH 2
+       -6.78E+00  4.21E+02
 
 TH 3
+        6.90E+00  1.98E+02  4.62E+02
 
 TH 4
+       -1.54E+01  3.16E+02 -3.04E+02  9.41E+02
 
 TH 5
+       -3.76E+00 -2.51E+02 -4.19E+02  2.73E+02  6.55E+02
 
 TH 6
+       -5.99E-01 -1.01E+00  1.89E+00 -4.16E+00 -1.01E+00  1.88E+02
 
 TH 7
+        3.99E-01  8.91E+00 -1.35E+01 -1.36E+01 -1.27E+01 -3.26E-01  1.27E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.41E+00 -2.02E+01 -3.19E+01  5.33E+01 -2.01E+00 -4.74E-01  1.68E+01  0.00E+00  4.49E+01
 
 TH10
+       -1.40E-01 -1.60E+01 -4.01E+01 -6.45E+00 -6.33E+01  1.60E-01  1.51E+01  0.00E+00  6.70E+00  8.18E+01
 
 TH11
+       -7.59E+00 -1.77E+01 -2.73E+01  1.47E-01 -2.98E+00  2.84E+00  9.62E+00  0.00E+00  4.38E+00  1.76E+01  1.75E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.847
Stop Time:
Wed Sep 29 16:13:48 CDT 2021
