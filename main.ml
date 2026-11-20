(* OPTION PRICER *)

type option_type= Call| Put
type market_region= European| American

let pay_off option_type k s=
  (*option_type: call or put; k=strike; s=spot*)
  match option_type with
  |Call -> max 0.0 (s -. k)
  |Put -> max 0.0 (k -. s)
;;

let arbre_option_pricing t s k r sigma q n option_type market_region=
  let p=Array.make (n+1) 0.0 in
  let deltaT =t/.(float_of_int n) in
  let up =exp(sigma*.sqrt(deltaT)) in
  let p0 =(up*.exp(-.q*.deltaT)-.exp(-.r*.deltaT))/.(up**2.-.1.) in
  let p1 =exp(-.r*.deltaT)-.p0 in 
  (* initial values at time T *)
  for i=0 to n do
    let spot=s*.up**(2.*.(float_of_int i)-.(float_of_int n)) in
    let value = pay_off option_type k spot in
    p.(i)<- value
  done;
  (* Move to earlier times *)
  for j=n-1 downto 0 do 
    for i=0 to j do
      let binomial_value=p0*.p.(i+1)+.p1*.p.(i) in 
      let spot=s*.up**(2.*.(float_of_int i)-.(float_of_int j)) in
      let exercise=pay_off option_type k spot in
      match market_region with 
      |American -> if binomial_value<exercise then p.(i)<-exercise else p.(i)<-binomial_value
      |European -> p.(i)<-binomial_value
    done;
  done;
  p.(0)
;;

let gaussian_generator ()=
  let u1= Random.float 1.0 in 
  let u2= Random.float 1.0 in 
  let u1 = if u1 = 0.0 then 1e-12 else u1 in
  (sqrt(-.2.*.log(u1))*.cos(2.*.Float.pi*.u2), sqrt(-.2.*.log(u1))*.sin(2.*.Float.pi*.u2))
;;


let monte_carlo_pricing opt_type t s k r sigma q m=
  let p=Array.make (2*m) 0.0 in
  for j=0 to m-1 do 
    let (z1, z2)= gaussian_generator () in
    let s1p=s*.exp((r-.q-.0.5*.sigma*.sigma)*.t+.(sigma*.sqrt(t)*.z1))
    and s1m=s*.exp((r-.q-.0.5*.sigma*.sigma)*.t+.(sigma*.sqrt(t)*.(-.z1))) 
    and s2p=s*.exp((r-.q-.0.5*.sigma*.sigma)*.t+.(sigma*.sqrt(t)*.z2))
    and s2m=s*.exp((r-.q-.0.5*.sigma*.sigma)*.t+.(sigma*.sqrt(t)*.(-.z2))) in
    let val1=((pay_off opt_type k s1p)+.(pay_off opt_type k s1m))/.2. 
    and val2=((pay_off opt_type k s2p)+.(pay_off opt_type k s2m))/.2. in
    p.(2*j)<-exp(-.r*.t)*.val1;
    p.(2*j+1)<-exp(-.r*.t)*.val2
  done;
  let moyenne= ref 0. in
  for j=0 to 2*m-1 do
    moyenne:=!moyenne+.p.(j)
  done;
  moyenne:=(1./.(float_of_int (2*m)))*.(!moyenne);
  let variance= ref 0. in
  for j=0 to 2*m-1 do
    variance:=(p.(j)-.(!moyenne))*.(p.(j)-.(!moyenne))+.(!variance)
  done;
  variance:=(1./.(float_of_int (2*m-1)))*.(!variance);
  (!moyenne, sqrt((!variance)/.(float_of_int (2*m))))
;;


(* Test with Black-Scholes *)

let erf_appox x=
  let sign=if x<0. then -1. else 1. in
  let x_abs=abs_float x in
  let p= 0.3275911 in 
  let a_1 = 0.254829592 in
  let a_2 = -0.284496736 in
  let a_3 = 1.421413741 in 
  let a_4 = -1.453152027 in
  let a_5 = 1.061405429 in
  let t=1./.(1.+.(p*.x_abs))in
  let poly = ((((a_5*.t+.a_4)*.t+.a_3)*.t+.a_2)*.t+.a_1)*.t in 
  sign*.(1.-.poly*.exp(-.x_abs*.x_abs))
;;

let phi x=
  let e= erf_appox (x /. sqrt 2.0) in 
  0.5*.(1.+.e)

  
;;

let black_scholes opt_type t s k r sigma q=
  let sqrt_t = sqrt t in
  let vol_sqrt_t = sigma *. sqrt_t in 
  let d_1 = (log(s/.k)+.(r-.q+.(0.5*.sigma*.sigma))*.t)/.vol_sqrt_t in
  let d_2= d_1-.vol_sqrt_t in 
  match opt_type with
  |Call-> s*.exp(-.q*.t)*.(phi d_1)-.k*.exp(-.r*.t)*.(phi d_2)
  |Put -> k*.exp(-.r*.t)*.(phi (-.d_2))-.s*.exp(-.q*.t)*.(phi (-.d_1))
;;

let t_, s_, k_, r_, simga_, q_, market_region_ = 0.5, 62.5, 62.5, 0.03, 0.048, 0.24, European
let test_black_call= black_scholes Call t_ s_ k_ r_ simga_ q_ 
let test_black_put= black_scholes Put t_ s_ k_ r_ simga_ q_ 
let (test_monte_carlo_call, err_call)= monte_carlo_pricing Call t_ s_ k_ r_ simga_ q_ 100000
let (test_monte_carlo_put, err_put)= monte_carlo_pricing Put t_ s_ k_ r_ simga_ q_  100000
let () =
  Printf.printf "Prix de l'option black Call: %.4f\n" test_black_call;
  Printf.printf "Prix de l'option black Put: %.4f\n" test_black_put;
  Printf.printf "----------------------------------------\n";
  Printf.printf "Prix de l'option monte carlo Call: %.4f\n" test_monte_carlo_call;
  Printf.printf "Prix de l'option monte carlo Put: %.4f\n" test_monte_carlo_put;
  Printf.printf "----------------------------------------\n";
  Printf.printf "err monte carlo Call: %.4f\n" (1.96*.err_call);
  Printf.printf "err monte carlo Put: %.4f\n" (1.96 *.err_put);
  Printf.printf "----------------------------------------\n";
  let l = [|10; 50; 200; 1000; 5000|] in
  for i = 0 to Array.length l - 1 do
    let n = l.(i) in
    let test_arbre_call = arbre_option_pricing t_ s_ k_ r_ simga_ q_ n Call market_region_ in
    let diff_call = abs_float (test_arbre_call -. test_black_call) in
    let test_arbre_put = arbre_option_pricing t_ s_ k_ r_ simga_ q_ n Put market_region_ in
    let diff_put = abs_float (test_arbre_put -. test_black_put) in
    Printf.printf "n = %d | Arbre Call: %.4f (Diff: %.4f) | Arbre Put: %.4f (Diff: %.4f)\n"
      n test_arbre_call diff_call test_arbre_put diff_put
  done


