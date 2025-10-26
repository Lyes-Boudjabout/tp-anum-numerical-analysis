funcprot(0);
clc;
clear;


mprintf("===========================================================\n");
mprintf("  TP ANUM - LYES BOUDJABOUT - HICHAM EDDHANE - G3 \n");
mprintf("  Étude numérique de la fonction : %s\n", "f(x) = exp(x) - x - 2");
mprintf("  Méthodes : Newton, Sécante, Halley\n");
mprintf("===========================================================\n");


// définition de la fonction
function y = f(x)
    y = exp(x) - x - 2;
endfunction

function y = df(x)
    y = exp(x) - 1;
endfunction

function y = d2f(x)
    y = exp(x);
endfunction

function y = g(x)
    y = cos(x) - x;
endfunction

function y = dg(x)
    y = - sin(x) - 1;
endfunction

function y = d2g(x)
    y = - cos(x);
endfunction


scf(0); // sélectionner la figure 0
clf();


 x = 0:0.001:5;
plot2d(x, f(x), 5); // style = 5
xgrid();
plot([0 5], [0 0], 2); // axe horizontal y=0
title("f(x) = exp(x) - x - 2");
xlabel("x");
ylabel("f(x)");


/************ Méthode de Newton ************/

/*** Fonction de Newton ***/
function [y,n_iter] = NEWTON(f, df, x0, e)
    n_iter = 0;
    r = x0;
    r2 = r + 0.01;
    
    while (abs(r - r2) > e)
        r2 = r;
        r = r2 - f(r2) / df(r2);
        n_iter = n_iter + 1;
    end
    y = r;
    mprintf("\n- Newton :\n");
    mprintf("  Valeur approchée de r : %f\n", y);
    mprintf("  Nombre d''itérations : %d\n", n_iter);
end

/*** Fin de Fonction de Newton ***/


/*** Test Newton **/
x0 = 1.3
e = 0.000001

[y,n_iter] = NEWTON(f, df, x0, e)

/*** Fin Test Newton **/

/************ Fin de Méthode de Newton ************/


/************ Méthode de La Sécante ************/

/*** Fonction de La Sécante ***/
function [y] = SECANTE(r0, r1, f, n)
 if f(r0)*f(r1) > 0 then
     mprintf("f(r0) et f(r1) sont de meme signe, on peut pas conrinuer\n")
     return
 end
 
 n_iter = 0
 e= 0.000001
 
 while((abs( f(r1)-f(r0) )>e) & (n_iter<n))
     temp = r1
     r1 = r1 - f(r1)*((r1-r0)/(f(r1)-f(r0)))
     r0 = temp
     n_iter = n_iter +1  
 end
 
 y = r1
 mprintf("\n- Sécante : \n")
 mprintf("la valeur approximative de r : %f\n",y)
 mprintf("nombre d''itérations : %d \n ",n_iter)
end
/*** Fin de Fonction de La Sécante ***/
/*** Test Sécante **/

r0 = 1;
r1 = 2;
n = 100;

[y] = SECANTE(r0, r1, f, n);
/*** Fin Test Sécante **/
/************ Fin de Méthode de La Sécante ************/



/************ Méthode de Halley ************/

/*** Fonction de Halley ***/
function [y,n_iter] = HALLEY(f, df, d2f, x0, e)

 n_iter = 0
 r = x0
 r2 = r+0.01
 
  while((abs(r-r2))>e) // critère d'arret = |rn+1-rn|<=e
   r2 = r // r2 == rn ; r == rn+1
   r = r2 - ( (2*f(r2)*df(r2)) / (2*df(r2)*df(r2)-f(r2)*d2f(r2)) ) // rn-[f(rn)/df(rn)] ==  g(rn) == rn+1
   n_iter=n_iter+1
  end

 y=r
 mprintf("\n- Halley : \n")
 mprintf("la valeur approximative de r : %f\n",y)
 mprintf("nombre d''itérations : %d\n\n\n",n_iter)

end
/*** Fin de Fonction de Halley ***/
/*** Test Halley **/

//x0 = 1.3;
e = 0.000001;

[y, n_iter] = HALLEY(f, df, d2f, x0, e);

/*** Fin Test Halley **/
/************ Fin de Méthode de Halley ************/

mprintf("==============================================================================================================\n");
mprintf("Commentaire sur les différences observées entre les méthodes et analyse de leur efficacité et avantages  \n");
mprintf("==============================================================================================================\n");

mprintf("\nObservation aprés l''éxecution : les 3 méthodes convergent vers la meme racine approximative\n")

mprintf("\n Newton :\n");
mprintf("  - Convergence quadratique {min ordre = 2} (très rapide si le point initial x0 est bien choisi).\n");
mprintf("  - Nécessite la dérivée première \n");
mprintf("  - Avantage : peu d’itérations lorsque x0 est proche de la racine.\n");
mprintf("  - Inconvénient : peut diverger si x0 est choisi loin de la racine.\n");


mprintf("\n Sécante :\n");
mprintf("  - Ne nécessite pas de dérivée \n");
mprintf("  - Convergence plus lente (ordre ≈ nombre or).\n");
mprintf("  - Avantage : utile quand la dérivée est difficile ou coûteuse à calculer.\n");
mprintf("  - Inconvénient : nécessite deux valeurs initiales bien choisies.\n");


mprintf("\n Halley :\n");
mprintf("  - Convergence au minimum quadratique.\n");
mprintf("  - Nécessite la 1ère et la 2ème dérivée (plus coûteux à évaluer).\n");
mprintf("  - Avantage : atteint la racine en très peu d’itérations.\n");
mprintf("  - Inconvénient : complexité plus élevée et 2 dérivées nécessaires.\n"); 


scf(1); // <- crée une nouvelle figure séparée
clf();

//x0 = 1.3;
x_vals = [x0];
errors = [];
for i = 1:10
    x_next = x_vals($) - f(x_vals($))/df(x_vals($));
    x_vals($+1) = x_next;
    errors($+1) = abs(x_next - x_vals($)); // Error between consecutive iterations
end

plot(0:length(x_vals)-1, x_vals, '-o', 'LineWidth', 2);
xgrid();
title("Méthode de Newton - Convergence", "fontsize", 4);
xlabel("Numéro d''itération (n)");
ylabel("x_n");

for i = 1:length(x_vals)
    xstring(i-1, x_vals(i) + 0.02, msprintf("%.4f", x_vals(i)));
end
xgrid();
title("Newton Method - Iteration Values", "fontsize", 3);
xlabel("Iteration Number (n)");
ylabel("x_n");
// Add value labels on points
for i = 1:length(x_vals)
    xstring(i-1, x_vals(i) + 0.02, msprintf("%.4f", x_vals(i)));
end


// Simple Secant Method Convergence Plot
scf(2); 
clf();

//r0 = 1;
//r1 = 2;
x_vals = [r0, r1];

// Calculate iterations
for i = 1:8
    r_next = r1 - f(r1)*((r1 - r0)/(f(r1) - f(r0)));
    x_vals($+1) = r_next;
    r0 = r1;
    r1 = r_next;
end

// Simple convergence plot
plot(0:length(x_vals)-1, x_vals, 's-', 'LineWidth', 2);
xgrid();
title("Secant Method - Iteration Values", "fontsize", 4);
xlabel("Iteration Number (n)");
ylabel("x_n");

// Add value labels
for i = 1:length(x_vals)
    if i <= 3 || i == length(x_vals)
        xstring(i-1, x_vals(i) + 0.02, msprintf("%.4f", x_vals(i)));
    end
end

// mprintf("Secant method converged to: %.6f\n", x_vals($));


// Simple Halley Method Convergence Plot - GUARANTEED TO WORK


// Calculate Halley iterations
//x0 = 1.3;
x_vals = [x0];

for i = 1:6
    x_current = x_vals($);
    numerator = 2 * f(x_current) * df(x_current);
    denominator = 2 * df(x_current)^2 - f(x_current) * d2f(x_current);
    x_next = x_current - numerator / denominator;
    x_vals($+1) = x_next;
end

// Simple plot - make it basic
scf(3); // Ensure we use figure 0
clf();

// Plot 1: Just the convergence
plot(0:length(x_vals)-1, x_vals, 'ro-');
xgrid();
title("Halley Method Convergence", "fontsize", 4);
xlabel("Iteration");
ylabel("x_n");

// Display the values
/* mprintf("\n=== Halley Method ===\n");
for i = 1:length(x_vals)
    mprintf("Iteration %d: x = %.8f\n", i-1, x_vals(i));
end */

scf(4); // sélectionner la figure 4
clf();


 x = 0:0.001:5;
plot2d(x, g(x), 5); // style = 5
xgrid();
plot([0 5], [0 0], 2); // axe horizontal y=0
title("g(x) = cos(x) - x");
xlabel("x");
ylabel("g(x)");

scf(5);// <- crée une nouvelle figure séparée
clf();

//x0 = 1.3;
x_vals = [x0];
errors = [];
for i = 1:10
    x_next = x_vals($) - g(x_vals($))/dg(x_vals($));
    x_vals($+1) = x_next;
    errors($+1) = abs(x_next - x_vals($)); // Error between consecutive iterations
end

plot(0:length(x_vals)-1, x_vals, '-o', 'LineWidth', 2);
xgrid();
title("Méthode de Newton - Convergence", "fontsize", 4);
xlabel("Numéro d''itération (n)");
ylabel("x_n");

for i = 1:length(x_vals)
    xstring(i-1, x_vals(i) + 0.02, msprintf("%.4f", x_vals(i)));
end
xgrid();
title("Newton Method - Iteration Values", "fontsize", 3);
xlabel("Iteration Number (n)");
ylabel("x_n");
// Add value labels on points
for i = 1:length(x_vals)
    xstring(i-1, x_vals(i) + 0.02, msprintf("%.4f", x_vals(i)));
end


// Simple Secant Method Convergence Plot
scf(6); 
clf();

//r0 = 1;
//r1 = 2;
x_vals = [r0, r1];

// Calculate iterations
for i = 1:8
    r_next = r1 - g(r1)*((r1 - r0)/(g(r1) - g(r0)));
    x_vals($+1) = r_next;
    r0 = r1;
    r1 = r_next;
end

// Simple convergence plot
plot(0:length(x_vals)-1, x_vals, 's-', 'LineWidth', 2);
xgrid();
title("Secant Method - Iteration Values", "fontsize", 4);
xlabel("Iteration Number (n)");
ylabel("x_n");

// Add value labels
for i = 1:length(x_vals)
    if i <= 3 || i == length(x_vals)
        xstring(i-1, x_vals(i) + 0.02, msprintf("%.4f", x_vals(i)));
    end
end

// mprintf("Secant method converged to: %.6f\n", x_vals($));


// Simple Halley Method Convergence Plot - GUARANTEED TO WORK


// Calculate Halley iterations
//x0 = 1.3;
x_vals = [x0];

for i = 1:6
    x_current = x_vals($);
    numerator = 2 * g(x_current) * dg(x_current);
    denominator = 2 * dg(x_current)^2 - g(x_current) * d2g(x_current);
    x_next = x_current - numerator / denominator;
    x_vals($+1) = x_next;
end

// Simple plot - make it basic
scf(7); // Ensure we use figure 0
clf();

// Plot 1: Just the convergence
plot(0:length(x_vals)-1, x_vals, 'ro-');
xgrid();
title("Halley Method Convergence", "fontsize", 4);
xlabel("Iteration");
ylabel("x_n");

// Display the values
/* mprintf("\n=== Halley Method ===\n");
for i = 1:length(x_vals)
    mprintf("Iteration %d: x = %.8f\n", i-1, x_vals(i));
end */

