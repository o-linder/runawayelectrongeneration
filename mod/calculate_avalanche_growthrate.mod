	  7  3   k820309    j          14.0        +K[_                                                                                                           
       src/calculate_avalanche_growthrate.f CALCULATE_AVALANCHE_GROWTHRATE              DP E_C E_CEFF_OVER_EC E_D ITER_ACC ITER_MAX LN_LAMBDA_0 LN_LAMBDA_C NU_EE NUBAR_D NUBAR_S PI                                                    
                            @                              
                            @                              
       PI                      @                              
       LN_LAMBDA_0 LN_LAMBDA_C                      @                              
       NU_EE NUBAR_D NUBAR_S                                                   KIND %         @                                                  
       #E_CEFF_OVER_E_C%SUM    #E_CEFF_OVER_E_C%LOG 	   #E_CEFF_OVER_E_C%SQRT 
   #E_CEFF_OVER_E_C%ABS    #Z    #Z0    #NJ    #NSPECIES    #T    #B                                                     SUM                                             	     LOG                                             
     SQRT                                                  ABS          
                                                          p          5 O p            5 O p                                   
                                                          p          5 O p            5 O p                                   
                                                     
    p          5 O p            5 O p                                    
                                                       
                                      
                
                                      
      %         @                                                  
       #AVALANCHE_GROWTHRATE_CLASSIC%SQRT    #AVALANCHE_GROWTHRATE_CLASSIC%SUM    #Z    #Z0    #NJ    #NSPECIES    #T    #EPAR    #EPS                                                    SQRT                                                 SUM          
                                                          p          5 � p        r        5 � p        r                               
                                                          p          5 � p        r        5 � p        r                               
                                                     
    p          5 � p        r        5 � p        r                                
                                                       
  @                                   
                
                                      
                
  @                                   
      %         @                                                  
       #AVALANCHE_GROWTHRATE%SQRT    #AVALANCHE_GROWTHRATE%SUM    #Z    #Z0 !   #NJ "   #NSPECIES     #T #   #EPAR $   #B %                                                   SQRT                                                 SUM          
  @                                                       p          5 � p        r         5 � p        r                                
  @                               !                        p          5 � p        r         5 � p        r                                
  @                              "                    
    p          5 � p        r         5 � p        r                                 
  @                                                     
  @                              #     
                
  @                              $     
                
  @                              %     
      %         @  @                            &                   
       #P_STAR%ABS '   #P_STAR%MAX (   #P_STAR%SQRT )   #P_STAR%SUM *   #Z +   #Z0 -   #NJ .   #NSPECIES ,   #T /   #EPAR 0   #B 1                                              '     ABS                                            (     MAX                                            )     SQRT                                            *     SUM          
  @                               +                        p          5 � p        r ,       5 � p        r ,                              
  @                               -                        p          5 � p        r ,       5 � p        r ,                              
  @                              .                    
 	   p          5 � p        r ,       5 � p        r ,                               
  @                               ,                     
  @                              /     
                
                                 0     
                
  @                              1     
         �   L      fn#fn 4   �   m   b   uapp(CALCULATE_AVALANCHE_GROWTHRATE    Y  @   J   DOUBLE     �  @   J   ELECTRIC_FIELDS #   �  C   J  PHYSICAL_CONSTANTS #     X   J  COULOMB_LOGARITHMS &   t  V   J  COLLISION_FREQUENCIES    �  =       KIND+DOUBLE 0     �       E_CEFF_OVER_E_C+ELECTRIC_FIELDS 4   �  <      E_CEFF_OVER_E_C%SUM+ELECTRIC_FIELDS 4   +  <      E_CEFF_OVER_E_C%LOG+ELECTRIC_FIELDS 5   g  =      E_CEFF_OVER_E_C%SQRT+ELECTRIC_FIELDS 4   �  <      E_CEFF_OVER_E_C%ABS+ELECTRIC_FIELDS 2   �  �   a   E_CEFF_OVER_E_C%Z+ELECTRIC_FIELDS 3   �  �   a   E_CEFF_OVER_E_C%Z0+ELECTRIC_FIELDS 3   (  �   a   E_CEFF_OVER_E_C%NJ+ELECTRIC_FIELDS 9   �  @   a   E_CEFF_OVER_E_C%NSPECIES+ELECTRIC_FIELDS 2     @   a   E_CEFF_OVER_E_C%T+ELECTRIC_FIELDS 2   L  @   a   E_CEFF_OVER_E_C%B+ELECTRIC_FIELDS -   �  �       AVALANCHE_GROWTHRATE_CLASSIC 2   h  =      AVALANCHE_GROWTHRATE_CLASSIC%SQRT 1   �  <      AVALANCHE_GROWTHRATE_CLASSIC%SUM /   �  �   a   AVALANCHE_GROWTHRATE_CLASSIC%Z 0   �	  �   a   AVALANCHE_GROWTHRATE_CLASSIC%Z0 0   I
  �   a   AVALANCHE_GROWTHRATE_CLASSIC%NJ 6   �
  @   a   AVALANCHE_GROWTHRATE_CLASSIC%NSPECIES /   =  @   a   AVALANCHE_GROWTHRATE_CLASSIC%T 2   }  @   a   AVALANCHE_GROWTHRATE_CLASSIC%EPAR 1   �  @   a   AVALANCHE_GROWTHRATE_CLASSIC%EPS %   �  �       AVALANCHE_GROWTHRATE *   �  =      AVALANCHE_GROWTHRATE%SQRT )     <      AVALANCHE_GROWTHRATE%SUM '   @  �   a   AVALANCHE_GROWTHRATE%Z (   �  �   a   AVALANCHE_GROWTHRATE%Z0 (   �  �   a   AVALANCHE_GROWTHRATE%NJ .   \  @   a   AVALANCHE_GROWTHRATE%NSPECIES '   �  @   a   AVALANCHE_GROWTHRATE%T *   �  @   a   AVALANCHE_GROWTHRATE%EPAR '     @   a   AVALANCHE_GROWTHRATE%B    \  �       P_STAR    *  <      P_STAR%ABS    f  <      P_STAR%MAX    �  =      P_STAR%SQRT    �  <      P_STAR%SUM      �   a   P_STAR%Z    �  �   a   P_STAR%Z0    �  �   a   P_STAR%NJ     7  @   a   P_STAR%NSPECIES    w  @   a   P_STAR%T    �  @   a   P_STAR%EPAR    �  @   a   P_STAR%B 