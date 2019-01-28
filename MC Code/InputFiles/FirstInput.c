type eigenPhysicsPackage; 

//pop      100000; 
pop      10000;
active   500; 
inactive 300; 
//seed     -6574747;
XSdata   ce2;

collisionOperator { type perNuclideCollisionOpCE; 
                   //type perMaterialCollisionOpMG; 
                  } 

transportOperator { //type transportOperatorST;
                     type transportOperatorDT;     
                  } 

inactiveTally { //clerk1 { type keffActiveClerk; display yes; } 
              } 

activeTally  { //clerk1 { type macroClerk; display no; map { type energyMap; grid log; min 1.0E-11; max 20.0; N 300; } }
               //clerk2 { type macroClerk; display no; map { type spaceMap; grid lin; axis x; min -10.0; max 10.0; N 20;} } 
               //clerk3 {type collProbClerk; display yes; materials (uo2 water); }
               //clerk4 {type dancoffBellClerk; XSmat uo2; fuelMat (uo2); modMat (water);
		       //Etop 0.02;
		       //Elow 4.0E-6;
		    //  map { type energyMap; grid log; min 0.5E-3; max 0.7E-3; N 50;}  
		      //map { type matXsMap; grid log; min 0.2; max 6.0; N 50; mat uo2;}
	            //  }
	      }

geometry { 
    type basicCellCSG;
    boundary (1 1 2 2 0 0);

    surfaces
    {
      squareBound { id 1; type zSquareCylinder; origin ( 0.0  0.0  0.0); halfwidth (10.0 10.0 0.0);}  
      zcyl        { id 2; type zCylinder      ; origin ( 0.0  0.0  0.0); radius 0.3;}

      planeX      { id 3; type xPlane         ; x 0.0;} 
      planeY      { id 4; type yPlane         ; y 0.0;} 
    }


    cells
    {

      fuel {id 1; surfaces (-2); filltype mat; mat uo2;  }
      mod  {id 2; surfaces (2) ; filltype mat; mat water; }
  
      topRight { id 10; surfaces ( 3  4 -1); filltype uni; universe 20 ;} 
      botRigth { id 11; surfaces ( 3 -4 -1); filltype uni; universe 21 ;}
      topLeft  { id 12; surfaces (-3  4 -1); filltype uni; universe 22 ;} 
      botLeft  { id 13; surfaces (-3 -4 -1); filltype uni; universe 23 ;}
  
      out     { id 3; surfaces (1 ); filltype outside;         }
      inside  { id 4; surfaces (-1); filltype uni; universe 17;}
    }

    universes
    {

      root
      {
	  id 1;
	  type cellUniverse; 
	  origin (0.0 0.0 0.0);
	//  cells (10 11 12 13 3);
	  cells ( 4 3); 
      }

      uni20 {id 20; type cellUniverse; origin (0.5 0.5 0.0); cells (1 2); }
      uni21 {id 21; type cellUniverse; origin (0.5 -0.5 0.0); cells (1 2); }
      uni22 {id 22; type cellUniverse; origin (-0.5 0.5 0.0); cells (1 2); }
      uni23 {id 23; type cellUniverse; origin (-0.5 -0.5 0.0); cells (1 2); }  

      uni30 {id 30; type cellUniverse; origin (0.0 0.0 0.0); cells (1 2); }
 
      uni31 {id 31; type pinUniverse; radii (0.3 0.0); fills (uo2 water); } 
 
  
      latUni 
      {   id 17; 
	  type   latUniverse; 
	  origin (0.0 0.0 0.0); 
	  pitch  (1.0 1.0 1.0); 
	  shape  (20   20   0); 
	  padMat water; 
	  map    ( 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31   
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31
	           31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31         
		   31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31 31);    
      }
    }   
} 


nuclearData {

  handles { 
     // mg1 isotropicMG;
     // mg2 transMG;
     // mg3 P1MG; 
     // ce1 byNucNoMT;
      ce2 byNucMT; 
  }
  
  
materials { 
    aceLibrary /home/mak60/myACE/JEF311.aceXS; 
    numberOfGroups 69; 

    water { 
      temp      273; 
      1001.03c    5.028E-02;
      8016.03c    2.505E-02;
      5010.03c    2.0E-005;
      //xsFile ./InputFiles/CAS8_12_1;  
      xsFile ./InputFiles/WIMS69HS_12_1;
      
    } 

    uo2 {  
      temp       273; 
      8016.03c   4.571E-02; 
      92238.03c  2.180E-02; 
      92235.03c  1.042E-03;
      //xsFile ./InputFiles/CAS8_11_1;
      xsFile ./InputFiles/WIMS69HS_11_1;      
    }	 

} 
  
}
  
  
  
