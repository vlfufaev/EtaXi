options(digits = 20);

bisection<-function(a, b, f, eps1){
  while((b - a) > eps1){
    if((f((b + a)/2))*(f(b)) > 0){
      b = ((b + a)/2);
    }
    else{
      a = ((b + a)/2);
    }
  }
  ab <- list("left" = a, "right" = b);
  return(ab);#
}

albet<-function(KAPPA, RELAMBDA, IMLAMBDA, eps1){
  K3 = -KAPPA/3;
  phi1 <-function(x){z <- (-1)*sqrt(x^2/3 + KAPPA/3 - RELAMBDA/3/x)*(3*x^2 + KAPPA - (x^2/3 + KAPPA/3 - RELAMBDA/3/x)) - IMLAMBDA}
  phi2 <-function(x){z <- sqrt(x^2/3 + KAPPA/3 - RELAMBDA/3/x)*(3*x^2 + KAPPA - (x^2/3 + KAPPA/3 - RELAMBDA/3/x))- IMLAMBDA}
  
  rootsalbet<-list('al1' = 0, 'al2' = 0, 'be1' = 0, 'be2' = 0, 'alp1' = 0, 'alp2' = 0, 'bet1' = 0, 'bet2' = 0);
  s = 1;
  while(((-sqrt(K3) - s*sqrt(K3)/3)^3 + KAPPA*(-sqrt(K3) - s*sqrt(K3)/3)) > RELAMBDA){
    s = s + 1;
  }
  l1 = bisection(-sqrt(K3) - s*sqrt(K3)/3, -sqrt(K3) - (s - 1)*sqrt(K3)/3, function(x){z <- x^3 + KAPPA*x - RELAMBDA}, eps1);
  rootsalbet$al1 =l1$left; rootsalbet$al2 = l1$right;

  s = 2;
  while(((-sqrt(K3)/s)^3 + KAPPA*(-sqrt(K3)/s)) > RELAMBDA){
    s = s + 1;
  }
  l1 = bisection(-sqrt(K3)/(s - 1), -sqrt(K3)/s, function(x){z <- x^3 + KAPPA*x - RELAMBDA}, eps1); 
  rootsalbet$be1 =l1$left; rootsalbet$be2 = l1$right;

  s = 1;
  while(phi1(rootsalbet$al1 - s*sqrt(K3)/3) > IMLAMBDA){
    s = s + 1;
  }
  l1 = bisection(rootsalbet$al1 - s*sqrt(K3)/3, rootsalbet$al1 - (s - 1)*sqrt(K3)/3, phi1, eps1); 
  rootsalbet$alp1 =l1$left; rootsalbet$alp2 = l1$right;
  s = 2;
  while(phi2(rootsalbet$be2/s) > IMLAMBDA){
    s = s + 1;
  }
  l1 = bisection(rootsalbet$be2/(s - 1), rootsalbet$be2/s, phi2, eps1);
  rootsalbet$bet1 =l1$left; rootsalbet$bet2 = l1$right;
  return(rootsalbet);
}

etaxi<-function(KAPPA, RELAMBDA, IMLAMBDA, eps1, eps2){
  eps1 = eps1/4;
  rootsalbet = albet(KAPPA, RELAMBDA, IMLAMBDA, eps1);
  al1 = rootsalbet$al1; al2 = rootsalbet$al2; alp1 = rootsalbet$alp1; alp2 = rootsalbet$alp2; 
  be1 = rootsalbet$be2; be2 = rootsalbet$be2; bet1 = rootsalbet$bet1; bet2 = rootsalbet$bet2;
  K3 = -KAPPA/3;
  S11 = -10; S12 = 10;n = 100000;#S(\alpha(\la, \kappa), \alpha(a, \kappa));
  while((abs(S11 - S12) > eps2)){
    S11 = 0; S12 = 0; 
    tmp = sqrt(-sqrt((alp2^3 + alp2*KAPPA - RELAMBDA)/3/alp2)*(8*alp2^2/3 + 2*KAPPA/3 + RELAMBDA/3/alp2) - IMLAMBDA);
    xnp = alp2;
    y1 = sqrt((alp2^3  - RELAMBDA)/(3*alp2) - K3);
    dx = (al1 - alp2)/n
  
    S11 = S11 + tmp*(y1 - sqrt((alp1^3 + KAPPA*alp1 - RELAMBDA)/(3*alp1)));#(s_0, s_1)
    S11 = S11 + sqrt(-IMLAMBDA)*(-sqrt((al1^3 + KAPPA*al1 - RELAMBDA)/(3*al1)));#(s_n, s_{n+1})
    S12 = S12 + sqrt(-sqrt((al1^3 + al1*KAPPA - RELAMBDA)/3/al1)*(8*al1^2/3 + 2*KAPPA/3 + RELAMBDA/3/al1) - IMLAMBDA)*(- sqrt((al1^3 + KAPPA*al1 - RELAMBDA)/(3*al1)));#(s_n, s_{n+1})
    
    for (s  in 0:(n - 1)){
      xnp = alp2 + (s + 1)*dx;#xnp = xn + dx;
      y2 = sqrt((xnp^3  - RELAMBDA)/(3*xnp) - K3);
      dy = y2 - y1;
    
      S12 = S12 + tmp*dy;
      tmp = sqrt(-y2*(3*xnp^2 + KAPPA - y2^2) - IMLAMBDA);
      S11 = S11 + tmp*dy;
      y1 = y2;
    }

    n = 10*n;
  }
  S21 = -10; S22 = 10;n = 100000;  #S(\beta(a, \kappa), \beta(\la, \kappa)); 
  while((abs(S21 - S22) > eps2)){
    S21 = 0; S22 = 0;
    tmp = sqrt(sqrt((be2^3 + KAPPA*be2 - RELAMBDA)/(3*be2))*(8*be2^2/3 + 2*KAPPA/3 + RELAMBDA/3/be2) - IMLAMBDA);
    xnp = be2;
    y1 = sqrt((be2^3  - RELAMBDA)/(3*be2) - K3);
    dx = (bet1 - be2)/n;
  
    S21 = S21 - sqrt(-IMLAMBDA)*y1;#(s_0, s_1)
    S22 = S22 - tmp*y1;#(s_0, s_1)  
    S21 = S21 - sqrt(sqrt((bet1^3 + KAPPA*bet1 - RELAMBDA)/3/bet1)*(8*bet1^2/3 + 2*KAPPA/3 + RELAMBDA/3/bet1) - IMLAMBDA)*(sqrt((bet2^3 + KAPPA*bet2 - RELAMBDA)/(3*bet2)) - sqrt((bet1^3 + KAPPA*bet1 - RELAMBDA)/(3*bet1)));#(s_n, s_{n+1})
  
    for (s  in 0:(n - 1)){
      xnp = be2 + (s + 1)*dx;
      y2 = sqrt((xnp^3  - RELAMBDA)/(3*xnp) - K3);
      dy =  y2 - y1;

      S21 = S21 - tmp*dy;
      tmp = sqrt(y2*(3*xnp^2 + KAPPA - y2^2) - IMLAMBDA);
      S22 = S22 - tmp*dy;
      y1 = y2;
    }  
    n = 10*n;
  } 
  I2 = IMLAMBDA^2;
  S31 = -10; S32 = 10;n = 100000;#S(-1, \alpha(a, \kappa));
  while((abs(S31 - S32) > eps2)){
    S32 = 0; S31 = 0;
    tmp = sqrt(sqrt(((-1)^3 + KAPPA*(-1) - RELAMBDA)^2 + I2) + IMLAMBDA)/sqrt(2);
    xnp = -1;
    S31 = S31 - sqrt(sqrt((al1^3 + KAPPA*al1 - RELAMBDA)^2 + I2) + IMLAMBDA)/sqrt(2)*(al2 - al1);#(s_n, s_{n+1})
    dx = (al1 + 1)/n;
  
    for (s  in 0:(n-1)){
      xnp = -1 + (s+1)*dx;
      S31 = S31 - tmp*dx;
      tmp = sqrt(sqrt((xnp^3 +KAPPA*xnp-RELAMBDA)^2 + I2) + IMLAMBDA)/sqrt(2);
      S32 = S32 - tmp*dx;
    }
    n = 10*n;
  } 
  S41 = -10; S42 = 10;n = 100000;  #S(\alpha(a, \kappa), \sqrt(KAPPA/3))
  while((abs(S41 - S42) > eps2)){
    S41 = 0; S42 = 0;
    tmp = sqrt(sqrt((al2^3 + KAPPA*al2 - RELAMBDA)^2 + I2) + IMLAMBDA)/sqrt(2);
    xnp = al2;
    S42 = S42 + tmp*(al2 - al1);#(s_0, s_1)
    dx = (-sqrt(K3) - al2)/n;
  
    for (s  in 0:(n-1)){
      xnp = al2 + (s + 1)*dx;
      S41 = S41 + tmp*dx;
      tmp = sqrt(sqrt((xnp^3 + KAPPA*xnp - RELAMBDA)^2 + I2) + IMLAMBDA)/sqrt(2);
      S42 = S42 + tmp*dx;
    }
    n = 10*n;
  } 
  S51 = -10; S52 = 10;n = 100000;  #S(\sqrt(KAPPA/3), \beta(a));
  while((abs(S51 - S52) > eps2)){
    S51 = 0; S52 = 0;
    dx = (be1 + sqrt(K3))/n;
    tmp = sqrt(sqrt(((-sqrt(K3))^3 + KAPPA*(-sqrt(K3)) - RELAMBDA)^2 + I2) + IMLAMBDA)/sqrt(2);
    xnp = -sqrt(K3);
  
    S52 = S52 + sqrt(sqrt((be1^3 +KAPPA*be1-RELAMBDA)^2 + I2) + IMLAMBDA)/sqrt(2)*(be2 - be1);#(s_n, s_{n+1})
    for (s  in 0:(n - 1)){
      xnp = -sqrt(K3) + (s + 1)*dx;
      S52 = S52 + tmp*dx;
      tmp = sqrt(sqrt((xnp^3 + KAPPA*xnp - RELAMBDA)^2 + I2) + IMLAMBDA)/sqrt(2);
      S51 = S51 + tmp*dx;
    } 
    n = 10*n;
  } 
  eta = S11 + S41 + S51 + S21; Eta = S12 + S42 + S52 + S22;xi = S11 - S31;Xi = S12 - S32;
  Etaxi<-list('xi' = xi, 'Xi' = Xi, 'eta' = eta, 'Eta' = Eta);
  return(Etaxi);
}