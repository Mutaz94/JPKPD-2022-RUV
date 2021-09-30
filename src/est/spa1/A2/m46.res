Wed Sep 29 23:22:28 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat46.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1239.83263684281        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7290E+02 -3.2044E+01  3.7165E+01 -5.6769E+01  1.1137E+02  4.8076E+01 -1.3488E+01  2.5380E-02  1.6901E+01 -5.3817E+01
            -1.6240E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1738.29176345538        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0268E+00  1.0738E+00  9.8139E-01  1.1006E+00  9.5087E-01  8.9552E-01  9.6043E-01  9.3075E-01  8.5300E-01  9.6542E-01
             2.4416E+00
 PARAMETER:  1.2647E-01  1.7124E-01  8.1211E-02  1.9589E-01  4.9625E-02 -1.0355E-02  5.9627E-02  2.8237E-02 -5.8992E-02  6.4810E-02
             9.9265E-01
 GRADIENT:   2.2886E+02  4.5346E+01 -2.5965E+00  7.6513E+01 -7.7139E+00 -1.1970E+01  7.2400E-01  6.8433E+00  3.1640E+00  1.0331E+01
             7.1436E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1748.69403203411        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0084E+00  8.8797E-01  4.3084E-01  1.1593E+00  5.5761E-01  9.1271E-01  1.3947E+00  1.8706E-01  8.1276E-01  4.4818E-01
             2.3243E+00
 PARAMETER:  1.0834E-01 -1.8819E-02 -7.4202E-01  2.4782E-01 -4.8409E-01  8.6682E-03  4.3267E-01 -1.5763E+00 -1.0732E-01 -7.0256E-01
             9.4341E-01
 GRADIENT:   1.5969E+02  3.4242E+01 -5.3484E+01  1.6744E+02  1.0429E+02 -5.1296E+00  3.1283E+01  2.6624E-01 -1.2090E-01  6.9361E-02
             5.9739E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1760.65207389417        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      337
 NPARAMETR:  9.5701E-01  7.0098E-01  3.8613E-01  1.2040E+00  4.6483E-01  9.1887E-01  1.4729E+00  8.8125E-02  7.6245E-01  5.6005E-01
             2.0865E+00
 PARAMETER:  5.6059E-02 -2.5528E-01 -8.5159E-01  2.8567E-01 -6.6608E-01  1.5385E-02  4.8720E-01 -2.3290E+00 -1.7122E-01 -4.7974E-01
             8.3547E-01
 GRADIENT:  -3.1757E+01  1.3656E+01 -6.1352E+01  1.1029E+02  8.2710E+01 -5.9614E+00  1.6321E+01  6.4207E-02 -1.5806E+01  2.4391E+00
            -4.3050E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1769.42505787958        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  9.5880E-01  4.8086E-01  3.0055E-01  1.1997E+00  3.3743E-01  9.4750E-01  1.4458E+00  1.5214E-02  8.5520E-01  5.9702E-01
             2.0028E+00
 PARAMETER:  5.7932E-02 -6.3218E-01 -1.1021E+00  2.8206E-01 -9.8639E-01  4.6069E-02  4.6867E-01 -4.0855E+00 -5.6415E-02 -4.1581E-01
             7.9454E-01
 GRADIENT:  -2.1009E+01 -2.4220E+00 -1.6683E+00  1.9720E-01  8.3059E+00  5.2447E+00  5.2000E+00  3.6229E-04 -1.2039E+00  2.9301E+00
            -1.3664E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1775.56311512277        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      691
 NPARAMETR:  9.6296E-01  3.6189E-01  1.7053E-01  1.0681E+00  2.2347E-01  9.5764E-01  6.5546E-01  1.0000E-02  1.0206E+00  7.1489E-01
             1.8672E+00
 PARAMETER:  6.2252E-02 -9.1641E-01 -1.6688E+00  1.6587E-01 -1.3985E+00  5.6713E-02 -3.2242E-01 -7.3092E+00  1.2039E-01 -2.3562E-01
             7.2442E-01
 GRADIENT:  -1.2304E+01 -3.8231E+00 -4.6782E+00 -1.4994E+01  1.7518E+01  5.5678E+00  2.8308E+00  0.0000E+00 -2.6853E+00 -9.6247E+00
            -4.3266E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1776.94363098019        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      867
 NPARAMETR:  9.6861E-01  3.2942E-01  1.4958E-01  1.0336E+00  2.0324E-01  9.4316E-01  4.4383E-01  1.0000E-02  1.0709E+00  7.5566E-01
             1.8668E+00
 PARAMETER:  6.8108E-02 -1.0104E+00 -1.7999E+00  1.3301E-01 -1.4934E+00  4.1480E-02 -7.1231E-01 -8.2758E+00  1.6853E-01 -1.8017E-01
             7.2423E-01
 GRADIENT:   1.1602E+00 -5.0692E+00 -5.4213E+00 -3.3095E+00  1.0122E+01  2.0255E-01  1.3301E+00  0.0000E+00  3.8952E-02 -1.3853E+00
             1.3678E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1777.58401831821        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1043
 NPARAMETR:  9.6814E-01  3.4921E-01  1.5935E-01  1.0542E+00  2.1243E-01  9.4299E-01  9.7711E-02  1.0000E-02  1.0525E+00  7.5983E-01
             1.8827E+00
 PARAMETER:  6.7617E-02 -9.5208E-01 -1.7367E+00  1.5281E-01 -1.4491E+00  4.1302E-02 -2.2257E+00 -8.4412E+00  1.5115E-01 -1.7466E-01
             7.3272E-01
 GRADIENT:  -6.4692E-01  1.7815E+00  2.1102E+00 -3.9905E-02 -3.6653E+00  7.4744E-02  5.7667E-02  0.0000E+00  5.5856E-01 -2.2639E-01
             1.8568E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1777.62747543733        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1221
 NPARAMETR:  9.6836E-01  3.4588E-01  1.5763E-01  1.0511E+00  2.1093E-01  9.4263E-01  2.8535E-02  1.0000E-02  1.0538E+00  7.6136E-01
             1.8818E+00
 PARAMETER:  6.7846E-02 -9.6167E-01 -1.7475E+00  1.4980E-01 -1.4562E+00  4.0916E-02 -3.4566E+00 -8.9440E+00  1.5237E-01 -1.7265E-01
             7.3223E-01
 GRADIENT:  -5.6483E-02 -1.0908E-02 -8.2843E-02  1.6936E-01  3.1313E-02 -2.4787E-02  4.9868E-03  0.0000E+00  2.3673E-02  9.1010E-04
             6.5358E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1777.62920046210        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1383
 NPARAMETR:  9.6838E-01  3.4598E-01  1.5769E-01  1.0511E+00  2.1099E-01  9.4268E-01  1.0000E-02  1.0000E-02  1.0535E+00  7.6137E-01
             1.8819E+00
 PARAMETER:  6.7871E-02 -9.6138E-01 -1.7471E+00  1.4981E-01 -1.4560E+00  4.0971E-02 -5.0093E+00 -9.4967E+00  1.5216E-01 -1.7263E-01
             7.3227E-01
 GRADIENT:   3.6512E-04  1.1360E-03  5.7979E-03 -5.7706E-03 -9.2116E-03 -1.6277E-04  0.0000E+00  0.0000E+00  2.7858E-03  2.0919E-04
            -1.6211E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1383
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1153E-03 -7.7952E-05  2.0646E-04 -7.9082E-03 -4.1346E-05
 SE:             2.9435E-02  1.5120E-04  2.5120E-04  2.7254E-02  2.7613E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6978E-01  6.0616E-01  4.1114E-01  7.7169E-01  9.9881E-01

 ETASHRINKSD(%)  1.3877E+00  9.9493E+01  9.9158E+01  8.6971E+00  7.4923E+00
 ETASHRINKVR(%)  2.7562E+00  9.9997E+01  9.9993E+01  1.6638E+01  1.4423E+01
 EBVSHRINKSD(%)  1.3674E+00  9.9470E+01  9.9265E+01  7.1392E+00  7.9292E+00
 EBVSHRINKVR(%)  2.7161E+00  9.9997E+01  9.9995E+01  1.3769E+01  1.5230E+01
 RELATIVEINF(%)  9.7241E+01  5.2594E-04  4.3316E-04  3.4032E+01  4.3652E+00
 EPSSHRINKSD(%)  3.0601E+01
 EPSSHRINKVR(%)  5.1837E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1777.6292004620991     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -858.69066725742641     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.71
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1777.629       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.68E-01  3.46E-01  1.58E-01  1.05E+00  2.11E-01  9.43E-01  1.00E-02  1.00E-02  1.05E+00  7.61E-01  1.88E+00
 


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
+        1.30E+03
 
 TH 2
+       -1.04E+01  2.54E+03
 
 TH 3
+       -5.75E+01  2.20E+03  2.34E+04
 
 TH 4
+       -3.69E+00  3.66E+01 -1.32E+03  7.85E+02
 
 TH 5
+        1.04E+02 -6.58E+03 -2.38E+04 -2.87E+02  3.66E+04
 
 TH 6
+        4.06E+00 -2.49E+01  5.38E+01 -1.02E+01 -2.01E-01  2.11E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.02E+01 -3.29E+01  3.11E+02 -1.66E+01  3.28E+01 -3.11E+00  0.00E+00  0.00E+00  1.30E+02
 
 TH10
+       -4.93E-01  1.52E+01  5.79E+01  1.30E+01 -3.00E+01  3.50E+00  0.00E+00  0.00E+00  5.38E+00  2.48E+02
 
 TH11
+       -1.74E+01 -1.34E+01 -9.51E+01 -5.10E+00  5.08E+01  2.99E+00  0.00E+00  0.00E+00  9.52E+00  1.23E+01  1.23E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.286
Stop Time:
Wed Sep 29 23:22:56 CDT 2021
