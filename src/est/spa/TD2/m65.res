Sat Sep 25 13:41:32 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat65.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1690.90394972794        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.8833E+01 -5.6089E+01 -5.6384E+01 -1.6189E+00  1.1489E+02  4.4571E+01 -5.9991E+00  5.5407E+00 -6.9634E+00 -1.8915E+01
            -3.6986E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1699.36022445844        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      133
 NPARAMETR:  1.0458E+00  1.0277E+00  1.0435E+00  1.0050E+00  9.4785E-01  8.3814E-01  1.0138E+00  9.7365E-01  1.0414E+00  1.0339E+00
             1.0131E+00
 PARAMETER:  1.4481E-01  1.2730E-01  1.4257E-01  1.0497E-01  4.6440E-02 -7.6575E-02  1.1369E-01  7.3300E-02  1.4056E-01  1.3332E-01
             1.1298E-01
 GRADIENT:   4.8431E+01 -2.9655E-01 -3.1114E+00  3.6751E+00  2.1573E-01 -2.2620E+01  3.2235E-01  3.3169E+00  8.0959E-01 -1.7962E+00
             9.5012E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1699.85120591393        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      312
 NPARAMETR:  1.0430E+00  1.0424E+00  1.0295E+00  9.9573E-01  9.5306E-01  8.4370E-01  9.6470E-01  8.3637E-01  1.0686E+00  1.0890E+00
             1.0111E+00
 PARAMETER:  1.4214E-01  1.4151E-01  1.2908E-01  9.5725E-02  5.1922E-02 -6.9953E-02  6.4057E-02 -7.8687E-02  1.6634E-01  1.8527E-01
             1.1100E-01
 GRADIENT:   3.9239E+01  2.0172E-01 -1.0028E+00  5.0498E+00 -5.4141E-01 -1.9652E+01 -7.5605E-01  8.8676E-01  2.5540E+00  3.2374E+00
             4.3008E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.94194815992        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      488
 NPARAMETR:  1.0307E+00  1.1927E+00  8.3477E-01  8.9108E-01  9.3250E-01  8.8748E-01  9.5165E-01  5.7242E-01  1.1164E+00  1.0304E+00
             1.0096E+00
 PARAMETER:  1.3023E-01  2.7626E-01 -8.0600E-02 -1.5324E-02  3.0114E-02 -1.9368E-02  5.0438E-02 -4.5789E-01  2.1014E-01  1.2998E-01
             1.0959E-01
 GRADIENT:  -2.3315E+00  2.4052E+00 -7.3495E-01  3.6494E+00 -3.0635E+00  7.1987E-01  4.0906E-01  4.9583E-01  2.2570E-01  1.0038E+00
             2.0493E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1701.26959901216        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      664
 NPARAMETR:  1.0331E+00  1.4837E+00  6.5791E-01  6.9901E-01  9.9513E-01  8.8722E-01  8.4067E-01  2.7868E-01  1.3086E+00  1.0349E+00
             1.0124E+00
 PARAMETER:  1.3259E-01  4.9451E-01 -3.1869E-01 -2.5809E-01  9.5123E-02 -1.9664E-02 -7.3550E-02 -1.1777E+00  3.6894E-01  1.3434E-01
             1.1228E-01
 GRADIENT:   1.4398E+00  6.4870E+00  4.8451E-01  2.7032E+00 -4.6571E+00  2.0917E-01  1.4498E-01  1.2564E-01 -8.3616E-02 -1.0737E-01
             1.8358E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1701.31475184859        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  1.0326E+00  1.5397E+00  6.3786E-01  6.5863E-01  1.0227E+00  8.8690E-01  8.2037E-01  2.1898E-01  1.3650E+00  1.0523E+00
             1.0133E+00
 PARAMETER:  1.3212E-01  5.3160E-01 -3.4964E-01 -3.1759E-01  1.2246E-01 -2.0018E-02 -9.8003E-02 -1.4188E+00  4.1113E-01  1.5103E-01
             1.1324E-01
 GRADIENT:   1.9020E-01 -5.0490E-01 -7.8739E-02 -4.6080E-01  3.0533E-02  9.4577E-02 -2.8450E-02  6.5175E-02  1.5884E-02 -3.9737E-02
             6.6614E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1701.34135126202        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1017
 NPARAMETR:  1.0323E+00  1.5298E+00  6.3641E-01  6.6497E-01  1.0173E+00  8.8613E-01  8.2579E-01  6.6085E-02  1.3534E+00  1.0503E+00
             1.0127E+00
 PARAMETER:  1.3179E-01  5.2511E-01 -3.5192E-01 -3.0802E-01  1.1719E-01 -2.0891E-02 -9.1412E-02 -2.6168E+00  4.0262E-01  1.4911E-01
             1.1259E-01
 GRADIENT:  -7.7770E-01 -7.8354E-01 -2.9518E-01  1.8779E-03  1.1405E+00 -2.6887E-01  6.4409E-02  4.5888E-03 -1.2364E-01  7.6472E-02
            -1.6220E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1701.34469844855        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1192
 NPARAMETR:  1.0326E+00  1.5319E+00  6.3358E-01  6.6360E-01  1.0161E+00  8.8673E-01  8.2480E-01  1.0000E-02  1.3557E+00  1.0483E+00
             1.0131E+00
 PARAMETER:  1.3208E-01  5.2653E-01 -3.5636E-01 -3.1008E-01  1.1593E-01 -2.0212E-02 -9.2617E-02 -4.7669E+00  4.0434E-01  1.4713E-01
             1.1300E-01
 GRADIENT:  -1.8184E-02 -1.7848E-02 -1.2141E-02 -1.5361E-03  3.4060E-02 -5.5176E-03  1.4004E-02  0.0000E+00 -6.6399E-03 -4.1048E-03
            -1.0054E-03

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1701.34469993537        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  1.0326E+00  1.5320E+00  6.3360E-01  6.6358E-01  1.0161E+00  8.8674E-01  8.2469E-01  1.0000E-02  1.3558E+00  1.0483E+00
             1.0131E+00
 PARAMETER:  1.3208E-01  5.2656E-01 -3.5634E-01 -3.1011E-01  1.1593E-01 -2.0199E-02 -9.2752E-02 -4.7614E+00  4.0442E-01  1.4717E-01
             1.1301E-01
 GRADIENT:  -2.3449E-03  2.5483E-04 -5.9921E-04  4.5968E-04  4.1393E-03 -2.9174E-04  6.2802E-04  0.0000E+00 -1.9617E-03 -1.7149E-03
            -7.5761E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1249
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2112E-05 -3.0589E-02 -3.0867E-04  2.1170E-02 -3.1850E-02
 SE:             2.9798E-02  2.1676E-02  1.1986E-04  2.3904E-02  2.3595E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9968E-01  1.5819E-01  1.0018E-02  3.7582E-01  1.7707E-01

 ETASHRINKSD(%)  1.7288E-01  2.7383E+01  9.9598E+01  1.9920E+01  2.0952E+01
 ETASHRINKVR(%)  3.4547E-01  4.7267E+01  9.9998E+01  3.5871E+01  3.7515E+01
 EBVSHRINKSD(%)  5.5813E-01  2.6597E+01  9.9664E+01  2.1178E+01  1.9232E+01
 EBVSHRINKVR(%)  1.1131E+00  4.6120E+01  9.9999E+01  3.7870E+01  3.4765E+01
 RELATIVEINF(%)  9.8783E+01  3.1537E+00  1.5425E-04  4.0979E+00  1.1957E+01
 EPSSHRINKSD(%)  4.4583E+01
 EPSSHRINKVR(%)  6.9290E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1701.3446999353669     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.19387337162868     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.68
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1701.345       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.53E+00  6.34E-01  6.64E-01  1.02E+00  8.87E-01  8.25E-01  1.00E-02  1.36E+00  1.05E+00  1.01E+00
 


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
+        1.31E+03
 
 TH 2
+       -6.82E+00  3.99E+02
 
 TH 3
+        1.08E+01  1.59E+02  3.52E+02
 
 TH 4
+       -1.73E+01  3.40E+02 -2.18E+02  8.60E+02
 
 TH 5
+       -3.94E+00 -2.26E+02 -3.44E+02  2.22E+02  6.40E+02
 
 TH 6
+       -1.88E-01 -1.84E+00  2.74E+00 -4.74E+00  2.78E-01  2.53E+02
 
 TH 7
+        4.90E-01  7.53E+00 -4.11E+00 -1.06E+01 -1.38E+01  1.11E+00  1.07E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.11E+00 -2.20E+01 -2.73E+01  4.72E+01  1.16E-01 -5.06E-01  1.91E+01  0.00E+00  5.10E+01
 
 TH10
+        3.26E-01 -1.52E+01 -3.93E+01 -6.00E+00 -5.23E+01  3.59E+00  1.40E+01  0.00E+00  5.47E+00  8.52E+01
 
 TH11
+       -8.83E+00 -1.76E+01 -2.96E+01  1.91E+00  2.61E-01  1.45E+00  1.01E+01  0.00E+00  4.24E+00  1.56E+01  2.13E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.695
Stop Time:
Sat Sep 25 13:41:58 CDT 2021
