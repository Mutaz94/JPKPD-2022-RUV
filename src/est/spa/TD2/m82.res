Wed Sep 29 19:20:51 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat82.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m82.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1620.24421375117        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2834E+02 -2.8898E+01 -6.4197E+01  7.9728E+01  1.3645E+02  6.5462E+01 -1.4741E+01  3.7951E+00  8.8901E+00 -2.5210E+01
            -5.5990E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1632.30815270147        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  1.0014E+00  1.0583E+00  1.0511E+00  9.7544E-01  9.5530E-01  9.3883E-01  1.0526E+00  9.8830E-01  9.9184E-01  1.0380E+00
             1.1482E+00
 PARAMETER:  1.0143E-01  1.5664E-01  1.4985E-01  7.5137E-02  5.4272E-02  3.6874E-02  1.5129E-01  8.8227E-02  9.1809E-02  1.3725E-01
             2.3816E-01
 GRADIENT:   7.0444E+00 -9.6673E+00 -7.5828E+00  3.6777E+00 -7.9423E-01  2.9546E-01 -9.1118E+00  3.5721E+00  3.5187E+00 -2.5275E+00
             7.5454E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1633.74700883989        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.9815E-01  1.1359E+00  1.0944E+00  9.3129E-01  1.0177E+00  9.4233E-01  1.1565E+00  8.8784E-01  9.2730E-01  1.1240E+00
             1.1335E+00
 PARAMETER:  9.8146E-02  2.2739E-01  1.9021E-01  2.8813E-02  1.1750E-01  4.0604E-02  2.4540E-01 -1.8965E-02  2.4518E-02  2.1686E-01
             2.2532E-01
 GRADIENT:  -5.8014E-01  1.4944E+00 -3.2279E-01  3.0410E+00  4.5187E+00  1.7156E+00 -1.7440E+00  6.4678E-01 -3.3846E+00  5.3384E-02
             7.4750E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1634.22633107697        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      530
 NPARAMETR:  1.0002E+00  1.2949E+00  8.0667E-01  8.2358E-01  9.6435E-01  9.3963E-01  1.0705E+00  3.8827E-01  1.0099E+00  1.0579E+00
             1.1224E+00
 PARAMETER:  1.0015E-01  3.5843E-01 -1.1484E-01 -9.4099E-02  6.3702E-02  3.7735E-02  1.6814E-01 -8.4606E-01  1.0985E-01  1.5629E-01
             2.1543E-01
 GRADIENT:  -4.8401E-01  9.0277E+00  1.8506E-01  8.5016E+00 -4.5256E+00 -4.1291E-01 -6.6155E-01  2.6293E-01 -2.5800E-01  1.7854E-02
            -1.4556E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1634.37343790001        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.0006E+00  1.4208E+00  7.3944E-01  7.3766E-01  1.0032E+00  9.4135E-01  1.0007E+00  2.1115E-01  1.0856E+00  1.0784E+00
             1.1265E+00
 PARAMETER:  1.0064E-01  4.5123E-01 -2.0186E-01 -2.0427E-01  1.0318E-01  3.9555E-02  1.0069E-01 -1.4552E+00  1.8212E-01  1.7543E-01
             2.1915E-01
 GRADIENT:  -1.0875E-01  3.1133E+00 -2.7389E-01  2.8839E+00 -4.1371E-01  2.3831E-01 -1.6472E-01  7.8908E-02 -3.6812E-01  7.1559E-02
             1.5808E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1634.40259154309        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  1.0006E+00  1.4417E+00  7.2831E-01  7.2341E-01  1.0101E+00  9.4088E-01  9.9015E-01  7.2181E-02  1.1061E+00  1.0817E+00
             1.1262E+00
 PARAMETER:  1.0064E-01  4.6584E-01 -2.1703E-01 -2.2378E-01  1.1008E-01  3.9062E-02  9.0098E-02 -2.5286E+00  2.0085E-01  1.7852E-01
             2.1883E-01
 GRADIENT:  -1.6783E-01  2.1480E+00 -1.2444E-01  2.2823E+00  2.6455E-01  2.8951E-02  1.1599E-01  8.8795E-03  1.4894E-01 -1.7546E-01
            -2.0718E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1634.41106254791        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1066             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0012E+00  1.4414E+00  7.2658E-01  7.2101E-01  1.0105E+00  9.4103E-01  9.8876E-01  1.5904E-02  1.1067E+00  1.0828E+00
             1.1261E+00
 PARAMETER:  1.0120E-01  4.6563E-01 -2.1940E-01 -2.2710E-01  1.1043E-01  3.9215E-02  8.8700E-02 -4.0412E+00  2.0137E-01  1.7955E-01
             2.1872E-01
 GRADIENT:   3.4301E+02  2.5959E+02  1.5143E+00  6.6878E+01  7.3800E+00  3.4519E+01  4.5278E+00  1.0946E-03  6.1231E+00  1.4317E+00
             1.7694E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1634.41140267524        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1205
 NPARAMETR:  1.0012E+00  1.4412E+00  7.2700E-01  7.2138E-01  1.0100E+00  9.4103E-01  9.8903E-01  1.0000E-02  1.1063E+00  1.0826E+00
             1.1260E+00
 PARAMETER:  1.0120E-01  4.6550E-01 -2.1883E-01 -2.2659E-01  1.0990E-01  3.9216E-02  8.8970E-02 -5.2903E+00  2.0099E-01  1.7935E-01
             2.1865E-01
 GRADIENT:   1.3026E+00 -6.8818E-01  1.3016E-01 -2.5365E-01 -1.5986E-01  1.0230E-01 -1.2450E-03  0.0000E+00  1.5259E-02  1.2587E-02
            -4.0962E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1634.41140267524        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1227
 NPARAMETR:  1.0012E+00  1.4412E+00  7.2700E-01  7.2138E-01  1.0100E+00  9.4103E-01  9.8903E-01  1.0000E-02  1.1063E+00  1.0826E+00
             1.1260E+00
 PARAMETER:  1.0120E-01  4.6550E-01 -2.1883E-01 -2.2659E-01  1.0990E-01  3.9216E-02  8.8970E-02 -5.2903E+00  2.0099E-01  1.7935E-01
             2.1865E-01
 GRADIENT:   1.3026E+00 -6.8818E-01  1.3016E-01 -2.5365E-01 -1.5986E-01  1.0230E-01 -1.2450E-03  0.0000E+00  1.5259E-02  1.2587E-02
            -4.0962E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1227
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.7595E-04 -1.6600E-02 -2.8718E-04  1.0145E-02 -2.6147E-02
 SE:             2.9778E-02  2.3124E-02  1.0871E-04  2.1714E-02  2.3386E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8993E-01  4.7284E-01  8.2488E-03  6.4037E-01  2.6354E-01

 ETASHRINKSD(%)  2.3922E-01  2.2532E+01  9.9636E+01  2.7255E+01  2.1654E+01
 ETASHRINKVR(%)  4.7787E-01  3.9987E+01  9.9999E+01  4.7081E+01  3.8619E+01
 EBVSHRINKSD(%)  6.0141E-01  2.1857E+01  9.9682E+01  2.9302E+01  1.9431E+01
 EBVSHRINKVR(%)  1.1992E+00  3.8936E+01  9.9999E+01  5.0017E+01  3.5086E+01
 RELATIVEINF(%)  9.8514E+01  2.7361E+00  1.0116E-04  2.1455E+00  1.0451E+01
 EPSSHRINKSD(%)  4.3057E+01
 EPSSHRINKVR(%)  6.7575E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1634.4114026752393     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -899.26057611150111     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.67
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1634.411       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.44E+00  7.27E-01  7.21E-01  1.01E+00  9.41E-01  9.89E-01  1.00E-02  1.11E+00  1.08E+00  1.13E+00
 


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
+        1.24E+03
 
 TH 2
+       -7.34E+00  3.48E+02
 
 TH 3
+        1.42E+01  9.27E+01  2.36E+02
 
 TH 4
+       -2.05E+01  3.55E+02 -2.10E+02  8.86E+02
 
 TH 5
+       -5.17E+00 -1.67E+02 -2.97E+02  2.52E+02  6.04E+02
 
 TH 6
+       -1.18E-01 -1.34E+00  2.55E+00 -4.35E+00 -1.34E+00  2.20E+02
 
 TH 7
+        1.06E+00  1.70E+01  8.09E-01 -1.81E+01 -6.15E+00 -8.15E-01  8.19E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.16E+00 -1.48E+01 -2.76E+01  3.86E+01  7.42E+00 -6.47E-01  2.19E+01  0.00E+00  5.05E+01
 
 TH10
+       -8.84E-01 -1.16E+01 -2.93E+01 -2.71E+00 -5.18E+01  4.49E-01  6.83E+00  0.00E+00  8.14E+00  7.27E+01
 
 TH11
+       -8.72E+00 -1.56E+01 -2.61E+01  3.66E+00  3.39E+00  3.36E+00  7.86E+00  0.00E+00  8.19E+00  1.67E+01  1.70E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.492
Stop Time:
Wed Sep 29 19:21:14 CDT 2021
