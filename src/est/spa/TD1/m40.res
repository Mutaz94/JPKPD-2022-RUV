Sat Sep 18 14:02:11 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat40.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1707.53330163823        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9475E+01 -1.4141E+01 -2.6993E+01  1.7815E+01  2.2869E+01  5.7683E+01 -6.3551E+00  7.0542E+00 -2.5569E+00  5.9832E+00
             1.2467E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1712.90227791599        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0178E+00  1.0513E+00  1.0612E+00  9.5618E-01  1.0593E+00  8.0563E-01  1.0711E+00  9.5231E-01  1.0268E+00  9.7479E-01
             9.7517E-01
 PARAMETER:  1.1766E-01  1.4998E-01  1.5939E-01  5.5193E-02  1.5758E-01 -1.1613E-01  1.6864E-01  5.1139E-02  1.2641E-01  7.4467E-02
             7.4852E-02
 GRADIENT:   7.8543E+01 -1.7729E+01 -1.0053E+01 -1.3936E+00  3.5584E+01 -1.8363E+01 -2.2929E+00  8.2804E-01  2.7043E+00 -6.6628E+00
            -1.4540E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1714.22848406557        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0121E+00  1.0771E+00  8.3867E-01  9.3135E-01  9.5718E-01  8.1950E-01  1.2502E+00  5.4459E-01  9.9570E-01  9.1411E-01
             9.7008E-01
 PARAMETER:  1.1200E-01  1.7429E-01 -7.5936E-02  2.8877E-02  5.6238E-02 -9.9061E-02  3.2330E-01 -5.0772E-01  9.5694E-02  1.0193E-02
             6.9624E-02
 GRADIENT:   5.1495E+01 -8.9011E-01 -1.1188E+01  5.1102E+00  1.9209E+01 -1.2086E+01  1.1912E+01  1.4971E+00  6.0647E+00  1.3407E+00
            -4.1155E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1715.47659305922        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.0142E+00  1.2263E+00  7.6156E-01  8.3825E-01  9.7893E-01  8.5845E-01  1.0316E+00  4.0810E-01  1.0764E+00  9.2144E-01
             9.7175E-01
 PARAMETER:  1.1409E-01  3.0398E-01 -1.7238E-01 -7.6439E-02  7.8708E-02 -5.2630E-02  1.3114E-01 -7.9623E-01  1.7366E-01  1.8185E-02
             7.1338E-02
 GRADIENT:   1.6713E+00 -5.4448E+00  7.4758E-01 -3.0256E+00 -1.6816E+00  2.4909E+00 -1.3230E+00  4.6136E-01  1.1639E+00 -2.8359E-01
            -9.8169E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1716.03397006377        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      471
 NPARAMETR:  1.0139E+00  1.4836E+00  6.4130E-01  6.7878E-01  1.0632E+00  8.5333E-01  8.9592E-01  1.2743E-01  1.2486E+00  9.6662E-01
             9.7640E-01
 PARAMETER:  1.1379E-01  4.9449E-01 -3.4425E-01 -2.8746E-01  1.6133E-01 -5.8614E-02 -9.9037E-03 -1.9602E+00  3.2202E-01  6.6055E-02
             7.6121E-02
 GRADIENT:  -1.6554E+00  3.8315E+00 -1.0238E+00  4.8977E+00  1.9819E+00 -4.1548E-04 -3.9071E-01  4.8196E-02 -5.6390E-01 -3.4761E-01
             6.6051E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1716.16676045727        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      646
 NPARAMETR:  1.0144E+00  1.6498E+00  5.7164E-01  5.6739E-01  1.1322E+00  8.5304E-01  8.2574E-01  3.1574E-02  1.4252E+00  1.0114E+00
             9.7691E-01
 PARAMETER:  1.1429E-01  6.0065E-01 -4.5925E-01 -4.6670E-01  2.2412E-01 -5.8948E-02 -9.1479E-02 -3.3554E+00  4.5431E-01  1.1132E-01
             7.6639E-02
 GRADIENT:  -3.4340E-01  2.0988E+00  4.6896E-01  3.4017E-01 -1.5328E+00 -9.8861E-02  1.7305E-01  2.6514E-03 -4.3237E-02  1.7405E-01
             2.9571E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1716.16998133722        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      826             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0145E+00  1.6667E+00  5.6285E-01  5.5534E-01  1.1404E+00  8.5326E-01  8.1868E-01  1.8158E-02  1.4466E+00  1.0156E+00
             9.7648E-01
 PARAMETER:  1.1442E-01  6.1083E-01 -4.7474E-01 -4.8817E-01  2.3135E-01 -5.8688E-02 -1.0007E-01 -3.9086E+00  4.6924E-01  1.1544E-01
             7.6194E-02
 GRADIENT:   5.5476E+01  7.0419E+01  3.7811E-01  1.2115E+01  1.4114E+00  4.0174E+00  7.1924E-01  9.6122E-04  1.9709E+00  1.6903E-01
             1.4161E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1716.17031674141        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  1.0145E+00  1.6667E+00  5.6279E-01  5.5535E-01  1.1404E+00  8.5326E-01  8.1854E-01  1.0000E-02  1.4466E+00  1.0154E+00
             9.7642E-01
 PARAMETER:  1.1441E-01  6.1085E-01 -4.7484E-01 -4.8816E-01  2.3141E-01 -5.8690E-02 -1.0023E-01 -4.6150E+00  4.6923E-01  1.1525E-01
             7.6141E-02
 GRADIENT:  -1.4741E-03 -2.7027E-02  3.9185E-03 -2.0927E-02 -5.4739E-02  9.0290E-04  2.6898E-02  0.0000E+00  5.7758E-03  3.1577E-02
             3.1236E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1716.17032137101        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1060
 NPARAMETR:  1.0145E+00  1.6667E+00  5.6279E-01  5.5535E-01  1.1405E+00  8.5326E-01  8.1845E-01  1.0000E-02  1.4466E+00  1.0153E+00
             9.7637E-01
 PARAMETER:  1.1441E-01  6.1086E-01 -4.7485E-01 -4.8815E-01  2.3142E-01 -5.8691E-02 -1.0034E-01 -4.6166E+00  4.6923E-01  1.1515E-01
             7.6085E-02
 GRADIENT:   8.4410E-04 -2.2629E-04  5.2458E-03  3.9759E-04 -2.2294E-02  2.4285E-04  6.2302E-03  0.0000E+00  7.7145E-05  7.8360E-03
             1.8230E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1060
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.1933E-05 -3.0924E-02 -2.9224E-04  2.6561E-02 -3.6800E-02
 SE:             2.9805E-02  2.3429E-02  1.0306E-04  2.2943E-02  2.2646E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9888E-01  1.8687E-01  4.5726E-03  2.4697E-01  1.0417E-01

 ETASHRINKSD(%)  1.4970E-01  2.1509E+01  9.9655E+01  2.3139E+01  2.4132E+01
 ETASHRINKVR(%)  2.9919E-01  3.8391E+01  9.9999E+01  4.0924E+01  4.2440E+01
 EBVSHRINKSD(%)  5.5563E-01  2.0703E+01  9.9701E+01  2.5023E+01  2.2272E+01
 EBVSHRINKVR(%)  1.1082E+00  3.7120E+01  9.9999E+01  4.3784E+01  3.9583E+01
 RELATIVEINF(%)  9.8793E+01  4.2319E+00  1.0831E-04  3.8285E+00  1.3125E+01
 EPSSHRINKSD(%)  4.4536E+01
 EPSSHRINKVR(%)  6.9237E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1716.1703213710102     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -981.01949480727205     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.56
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1716.170       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.67E+00  5.63E-01  5.55E-01  1.14E+00  8.53E-01  8.18E-01  1.00E-02  1.45E+00  1.02E+00  9.76E-01
 


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
+        1.47E+03
 
 TH 2
+       -7.13E+00  3.97E+02
 
 TH 3
+        1.20E+01  1.50E+02  3.53E+02
 
 TH 4
+       -1.82E+01  3.34E+02 -2.79E+02  9.84E+02
 
 TH 5
+       -5.69E+00 -1.77E+02 -2.96E+02  2.53E+02  5.02E+02
 
 TH 6
+       -1.62E+00 -1.67E+00  7.94E-01 -2.18E+00 -2.92E+00  2.64E+02
 
 TH 7
+       -4.52E+00  8.75E+00 -1.14E+01 -1.42E+01 -1.33E+01  2.66E+00  1.16E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.20E+00 -1.93E+01 -3.10E+01  5.63E+01  3.50E-01 -4.23E-01  1.89E+01  0.00E+00  3.96E+01
 
 TH10
+       -8.02E-01 -1.45E+01 -3.60E+01 -3.45E+00 -5.79E+01 -2.52E+00  1.13E+01  0.00E+00  5.02E+00  7.93E+01
 
 TH11
+       -1.15E+01 -1.78E+01 -2.78E+01  2.82E+00  2.42E+00  3.38E+00  1.09E+01  0.00E+00  4.04E+00  1.74E+01  2.22E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.750
Stop Time:
Sat Sep 18 14:02:31 CDT 2021
