Fri Sep 24 19:57:11 CDT 2021
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
$DATA ../../../../data/int/B/dat100.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3812.03581764003        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3557E+02 -5.3439E+01 -1.0168E+01  3.7215E+01  1.8952E+01 -1.3475E+01  2.2302E+01 -6.6082E+00  1.4272E+01  3.7772E+00
            -9.8508E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3821.32632767310        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7714E-01  1.2181E+00  1.1745E+00  8.8267E-01  1.2225E+00  1.0637E+00  7.1042E-01  1.1142E+00  8.8195E-01  1.0980E+00
             1.0275E+00
 PARAMETER:  7.6874E-02  2.9733E-01  2.6085E-01 -2.4799E-02  3.0093E-01  1.6179E-01 -2.4190E-01  2.0817E-01 -2.5617E-02  1.9349E-01
             1.2713E-01
 GRADIENT:   7.6788E+01  2.1583E+01  1.2566E+01 -5.9447E+00  6.2452E+01  1.7641E+01 -1.3797E+01 -1.9912E+01 -3.2403E+01 -1.9149E+01
            -5.1460E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3822.17455172895        NO. OF FUNC. EVALS.: 153
 CUMULATIVE NO. OF FUNC. EVALS.:      240
 NPARAMETR:  9.8171E-01  1.2236E+00  1.1692E+00  8.7999E-01  1.2264E+00  1.0650E+00  6.9461E-01  1.1432E+00  8.9766E-01  1.1108E+00
             1.0330E+00
 PARAMETER:  8.1540E-02  3.0177E-01  2.5633E-01 -2.7844E-02  3.0408E-01  1.6296E-01 -2.6440E-01  2.3380E-01 -7.9585E-03  2.0510E-01
             1.3249E-01
 GRADIENT:   4.3928E+01 -4.4987E+00  7.4138E+00 -1.1412E+01  4.9881E+01  9.2033E+00 -1.4855E+01 -1.7566E+01 -3.0471E+01 -1.8379E+01
            -3.7283E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3822.25668222254        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      420
 NPARAMETR:  9.8154E-01  1.2238E+00  1.1714E+00  8.8073E-01  1.2270E+00  1.0595E+00  6.9370E-01  1.1496E+00  8.9901E-01  1.1075E+00
             1.0329E+00
 PARAMETER:  8.1365E-02  3.0198E-01  2.5823E-01 -2.7006E-02  3.0456E-01  1.5779E-01 -2.6572E-01  2.3940E-01 -6.4662E-03  2.0207E-01
             1.3232E-01
 GRADIENT:   4.4020E+01 -3.1849E+00  7.5747E+00 -9.4584E+00  5.0365E+01  7.2585E+00 -1.4976E+01 -1.7397E+01 -3.0062E+01 -1.9268E+01
            -3.7348E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3822.29132852143        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  9.8154E-01  1.2234E+00  1.1714E+00  8.8176E-01  1.2270E+00  1.0553E+00  6.9370E-01  1.1496E+00  8.9901E-01  1.1075E+00
             1.0329E+00
 PARAMETER:  8.1365E-02  3.0164E-01  2.5823E-01 -2.5839E-02  3.0456E-01  1.5383E-01 -2.6572E-01  2.3940E-01 -6.4662E-03  2.0207E-01
             1.3232E-01
 GRADIENT:   4.4354E+01 -2.4795E+00  7.4065E+00 -7.4553E+00  5.0995E+01  5.7456E+00 -1.4955E+01 -1.7393E+01 -2.9987E+01 -1.9249E+01
            -3.7302E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3822.32373658653        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      777
 NPARAMETR:  9.8154E-01  1.2228E+00  1.1714E+00  8.8307E-01  1.2270E+00  1.0500E+00  6.9370E-01  1.1496E+00  8.9901E-01  1.1075E+00
             1.0329E+00
 PARAMETER:  8.1365E-02  3.0118E-01  2.5823E-01 -2.4349E-02  3.0456E-01  1.4877E-01 -2.6572E-01  2.3940E-01 -6.4662E-03  2.0207E-01
             1.3232E-01
 GRADIENT:   4.4784E+01 -1.6387E+00  7.1891E+00 -4.9276E+00  5.1820E+01  3.7950E+00 -1.4927E+01 -1.7388E+01 -2.9896E+01 -1.9225E+01
            -3.7240E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3823.27136667108        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      962            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8149E-01  1.2238E+00  1.1716E+00  8.8303E-01  1.1836E+00  1.0398E+00  6.9379E-01  1.1497E+00  8.9901E-01  1.1076E+00
             1.0329E+00
 PARAMETER:  8.1315E-02  3.0198E-01  2.5836E-01 -2.4399E-02  2.6856E-01  1.3908E-01 -2.6558E-01  2.3952E-01 -6.4662E-03  2.0217E-01
             1.3239E-01
 GRADIENT:   8.6951E+01  5.5531E+01  2.0528E+01 -8.5166E+00  1.1349E+01  7.2036E+00 -1.3292E+01 -1.6490E+01 -2.9910E+01 -1.4886E+01
            -3.5874E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3824.32302114898        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1143             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7979E-01  1.2150E+00  1.1716E+00  8.8385E-01  1.1808E+00  1.0370E+00  7.0299E-01  1.1685E+00  9.0777E-01  1.1124E+00
             1.0342E+00
 PARAMETER:  7.9584E-02  2.9477E-01  2.5835E-01 -2.3465E-02  2.6623E-01  1.3638E-01 -2.5242E-01  2.5574E-01  3.2331E-03  2.0653E-01
             1.3358E-01
 GRADIENT:   8.3596E+01  4.1176E+01  1.9044E+01 -1.4947E+01  1.3117E+01  6.1805E+00 -1.0709E+01 -1.4778E+01 -2.6220E+01 -1.2334E+01
            -3.0979E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3824.32395231353        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1307
 NPARAMETR:  9.7979E-01  1.2150E+00  1.1716E+00  8.8385E-01  1.1808E+00  1.0393E+00  7.0299E-01  1.1685E+00  9.0777E-01  1.1124E+00
             1.0342E+00
 PARAMETER:  7.9584E-02  2.9477E-01  2.5835E-01 -2.3465E-02  2.6623E-01  1.3854E-01 -2.5242E-01  2.5574E-01  3.2331E-03  2.0653E-01
             1.3358E-01
 GRADIENT:   4.2154E+01  1.3636E+01  1.7928E+01 -2.1992E+01  1.7759E+00 -1.9408E-04 -1.1802E+01 -1.4931E+01 -2.6841E+01 -1.3400E+01
            -3.1216E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1307
 NO. OF SIG. DIGITS IN FINAL EST.:  4.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.9319E-02 -3.6253E-02 -3.1836E-02  2.6936E-02 -2.4883E-02
 SE:             2.9832E-02  2.1304E-02  1.8922E-02  2.8752E-02  2.7145E-02
 N:                     100         100         100         100         100

 P VAL.:         5.1725E-01  8.8808E-02  9.2471E-02  3.4885E-01  3.5931E-01

 ETASHRINKSD(%)  5.9265E-02  2.8630E+01  3.6610E+01  3.6763E+00  9.0604E+00
 ETASHRINKVR(%)  1.1849E-01  4.9063E+01  5.9817E+01  7.2175E+00  1.7300E+01
 EBVSHRINKSD(%)  2.6472E-01  3.2634E+01  4.5410E+01  1.2502E+01  1.1593E+01
 EBVSHRINKVR(%)  5.2875E-01  5.4619E+01  7.0200E+01  2.3441E+01  2.1843E+01
 RELATIVEINF(%)  9.9470E+01  1.8040E+01  2.4101E+01  3.7066E+01  3.6861E+01
 EPSSHRINKSD(%)  1.9458E+01
 EPSSHRINKVR(%)  3.5130E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3824.3239523135298     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2170.2345925451191     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.08
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.81
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3824.324       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.22E+00  1.17E+00  8.84E-01  1.18E+00  1.04E+00  7.03E-01  1.17E+00  9.08E-01  1.11E+00  1.03E+00
 


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
+        3.30E+09
 
 TH 2
+       -2.71E+00  4.94E+08
 
 TH 3
+        2.85E-01  2.79E+01  6.91E+08
 
 TH 4
+       -8.48E-01  5.32E+02  1.18E+09  4.05E+09
 
 TH 5
+       -7.63E+02  2.81E+08 -1.20E+02  1.14E+09  3.20E+08
 
 TH 6
+       -4.23E+01 -8.38E-01  1.83E-01  9.51E-01 -7.21E+02  1.86E+02
 
 TH 7
+       -1.78E+00 -4.98E+08 -9.17E+00  3.19E+00 -5.68E+08 -2.33E+00  2.01E+09
 
 TH 8
+       -1.08E+09 -7.98E+00 -8.49E+01  3.28E+01 -1.69E+03 -2.75E-01  4.44E+00  3.54E+08
 
 TH 9
+        3.29E+00  9.24E+00  1.15E+09  3.95E+09  1.11E+09 -1.20E+00  1.97E+09  9.64E+01  3.84E+09
 
 TH10
+        1.41E+09 -2.43E+01  2.15E+00  1.17E+01 -4.38E+08  1.26E+00  2.98E+01 -4.61E+08 -2.85E-01  6.00E+08
 
 TH11
+       -1.74E+03 -4.62E+01 -2.09E+01  1.02E+01 -8.85E+00 -1.64E+03 -1.29E+09  6.36E+05  2.51E+01  5.42E+02  3.32E+09
 
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
 #CPUT: Total CPU Time in Seconds,       48.987
Stop Time:
Fri Sep 24 19:58:03 CDT 2021
