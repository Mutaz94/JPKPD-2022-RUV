Sat Sep 18 12:32:28 CDT 2021
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
 GRADIENT:   8.5922E+01 -3.7115E+01 -3.0605E+01 -1.4395E+01  5.9085E+01 -3.4763E+00 -2.9841E+00  4.7874E+00 -8.2808E+00 -2.9082E+00
            -4.5663E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1622.76884071533        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  9.8029E-01  1.0212E+00  1.0265E+00  1.0112E+00  9.6894E-01  1.0080E+00  1.0033E+00  9.8706E-01  1.0220E+00  9.9215E-01
             1.0907E+00
 PARAMETER:  8.0092E-02  1.2102E-01  1.2612E-01  1.1115E-01  6.8452E-02  1.0794E-01  1.0334E-01  8.6980E-02  1.2175E-01  9.2122E-02
             1.8685E-01
 GRADIENT:  -1.2889E+00  4.4379E+00 -2.7740E-01  3.2876E+00 -6.2443E+00 -2.6012E+00 -2.5696E-02  2.8974E+00 -3.4140E+00  2.8307E+00
            -5.2972E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1622.95071154039        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      313
 NPARAMETR:  9.7910E-01  1.0373E+00  1.0045E+00  1.0006E+00  9.7042E-01  1.0147E+00  9.7474E-01  8.5693E-01  1.0377E+00  9.9243E-01
             1.0913E+00
 PARAMETER:  7.8878E-02  1.3658E-01  1.0444E-01  1.0062E-01  6.9977E-02  1.1459E-01  7.4411E-02 -5.4398E-02  1.3702E-01  9.2400E-02
             1.8734E-01
 GRADIENT:  -4.2622E+00  5.1220E+00  3.1586E+00  4.6753E+00 -2.2279E+00 -3.1607E-02 -1.4573E+00 -7.3680E-01 -3.1893E+00  2.3179E-01
            -5.9597E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1623.26689514081        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  9.8354E-01  1.2273E+00  8.2747E-01  8.7169E-01  9.7541E-01  1.0166E+00  9.2054E-01  7.2534E-01  1.1413E+00  9.4963E-01
             1.1131E+00
 PARAMETER:  8.3405E-02  3.0478E-01 -8.9380E-02 -3.7317E-02  7.5098E-02  1.1642E-01  1.7206E-02 -2.2112E-01  2.3220E-01  4.8316E-02
             2.0713E-01
 GRADIENT:   1.6718E+00  4.5279E+00  1.7575E+00  3.5794E+00 -2.7811E+00  9.9511E-02  4.3234E-01 -4.7877E-02  1.3319E-01 -2.7575E-01
             2.8233E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1623.53195784101        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  9.8408E-01  1.4769E+00  6.2579E-01  7.0053E-01  1.0092E+00  1.0175E+00  8.3514E-01  4.9410E-01  1.2991E+00  9.3812E-01
             1.1035E+00
 PARAMETER:  8.3954E-02  4.8996E-01 -3.6873E-01 -2.5592E-01  1.0920E-01  1.1735E-01 -8.0161E-02 -6.0501E-01  3.6165E-01  3.6118E-02
             1.9845E-01
 GRADIENT:   5.5671E-01  1.3479E+00 -7.8283E-01  1.7302E+00  6.2579E-01 -7.4624E-02  2.8380E-01  2.7069E-01 -4.7351E-01 -6.2939E-03
            -5.1750E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1623.58858737368        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      841
 NPARAMETR:  9.8393E-01  1.5904E+00  5.5752E-01  6.2327E-01  1.0412E+00  1.0178E+00  7.9734E-01  3.6782E-01  1.4037E+00  9.5107E-01
             1.1056E+00
 PARAMETER:  8.3803E-02  5.6399E-01 -4.8425E-01 -3.7277E-01  1.4035E-01  1.1762E-01 -1.2647E-01 -9.0015E-01  4.3915E-01  4.9829E-02
             2.0035E-01
 GRADIENT:  -6.4408E-02  9.1637E-02 -5.0927E-01  1.9471E-01  3.4019E-01 -7.6851E-03 -1.5265E-01  2.0963E-01  2.8347E-02  1.3228E-01
             4.3628E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1623.65886635085        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1018
 NPARAMETR:  9.8407E-01  1.5913E+00  5.3493E-01  6.2060E-01  1.0290E+00  1.0179E+00  8.0151E-01  1.4886E-01  1.3972E+00  9.3597E-01
             1.1051E+00
 PARAMETER:  8.3945E-02  5.6452E-01 -5.2561E-01 -3.7706E-01  1.2855E-01  1.1774E-01 -1.2126E-01 -1.8047E+00  4.3443E-01  3.3829E-02
             1.9995E-01
 GRADIENT:   5.7434E-02 -1.1980E+00 -7.7990E-01  3.4724E-01  1.8070E+00 -1.1486E-02 -4.2578E-03  3.1199E-02 -1.8886E-01 -1.4026E-01
             1.8413E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1623.67388674683        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1193
 NPARAMETR:  9.8408E-01  1.5970E+00  5.2998E-01  6.1677E-01  1.0284E+00  1.0180E+00  8.0040E-01  2.6267E-02  1.4037E+00  9.3603E-01
             1.1048E+00
 PARAMETER:  8.3951E-02  5.6815E-01 -5.3492E-01 -3.8326E-01  1.2801E-01  1.1781E-01 -1.2264E-01 -3.5394E+00  4.3913E-01  3.3888E-02
             1.9966E-01
 GRADIENT:   6.5239E-02  1.6907E-01 -7.6571E-02  1.3315E-01 -3.0928E-02  9.2558E-03  6.9703E-02  9.2506E-04  6.6916E-03  1.8367E-02
             9.2461E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1623.67433984370        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1326
 NPARAMETR:  9.8406E-01  1.5977E+00  5.3018E-01  6.1610E-01  1.0292E+00  1.0179E+00  7.9975E-01  1.0000E-02  1.4050E+00  9.3671E-01
             1.1046E+00
 PARAMETER:  8.3912E-02  5.6867E-01 -5.3461E-01 -3.8426E-01  1.2874E-01  1.1779E-01 -1.2340E-01 -4.9293E+00  4.4007E-01  3.4719E-02
             1.9946E-01
 GRADIENT:  -1.9379E-02  7.7197E-02 -8.7607E-03  3.0287E-02 -1.5143E-02  1.6785E-03  4.3263E-03  0.0000E+00  3.6189E-04  7.2938E-03
            -7.5399E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1326
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0573E-04 -2.9587E-02 -2.8726E-04  2.2052E-02 -3.4080E-02
 SE:             2.9821E-02  2.2787E-02  1.1586E-04  2.3756E-02  2.2113E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9450E-01  1.9414E-01  1.3158E-02  3.5326E-01  1.2328E-01

 ETASHRINKSD(%)  9.5380E-02  2.3661E+01  9.9612E+01  2.0414E+01  2.5917E+01
 ETASHRINKVR(%)  1.9067E-01  4.1724E+01  9.9998E+01  3.6661E+01  4.5117E+01
 EBVSHRINKSD(%)  5.0463E-01  2.3189E+01  9.9668E+01  2.1406E+01  2.4912E+01
 EBVSHRINKVR(%)  1.0067E+00  4.1001E+01  9.9999E+01  3.8230E+01  4.3618E+01
 RELATIVEINF(%)  9.8959E+01  4.0313E+00  1.3858E-04  4.7187E+00  9.7374E+00
 EPSSHRINKSD(%)  4.3951E+01
 EPSSHRINKVR(%)  6.8585E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1623.6743398437004     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -888.52351327996223     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.36
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1623.674       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.84E-01  1.60E+00  5.30E-01  6.16E-01  1.03E+00  1.02E+00  8.00E-01  1.00E-02  1.41E+00  9.37E-01  1.10E+00
 


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
+       -6.49E+00  4.22E+02
 
 TH 3
+        6.46E+00  2.00E+02  4.69E+02
 
 TH 4
+       -1.48E+01  3.16E+02 -2.99E+02  9.28E+02
 
 TH 5
+       -4.11E+00 -2.57E+02 -4.28E+02  2.71E+02  6.72E+02
 
 TH 6
+       -3.90E-02 -9.62E-01  2.55E+00 -3.75E+00 -1.45E+00  1.89E+02
 
 TH 7
+        1.32E+00  8.92E+00 -1.25E+01 -1.49E+01 -1.19E+01 -7.74E-01  1.25E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.56E+00 -2.06E+01 -3.14E+01  5.20E+01 -2.15E+00 -2.40E-01  1.67E+01  0.00E+00  4.68E+01
 
 TH10
+       -1.99E-01 -1.60E+01 -4.04E+01 -8.73E+00 -6.29E+01  7.45E-01  1.69E+01  0.00E+00  7.19E+00  8.16E+01
 
 TH11
+       -8.93E+00 -1.77E+01 -2.78E+01 -3.12E-01 -2.94E+00  2.46E+00  9.60E+00  0.00E+00  4.45E+00  1.95E+01  1.75E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.348
Stop Time:
Sat Sep 18 12:32:52 CDT 2021
