// https://stackoverflow.com/questions/29001275/roots-of-a-polynomial-mod-a-prime/29022281#29022281
// https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields

// https://www.imsc.res.in/~vikram/poly-2012/euclid.pdf
// https://www.csd.uwo.ca/~mmorenom/CS874/Lectures/Newton2Hensel.html/node12.html
// https://ja.wikipedia.org/wiki/%E3%82%B7%E3%83%A7%E3%83%BC%E3%83%B3%E3%83%8F%E3%83%BC%E3%82%B2%E3%83%BB%E3%82%B9%E3%83%88%E3%83%A9%E3%83%83%E3%82%BB%E3%83%B3%E6%B3%95
// https://algo.inria.fr/seminars/sem08-09/kruppa-slides.pdf
// https://arxiv.org/abs/1602.04562

// https://ja.wikipedia.org/wiki/%E3%83%A2%E3%83%B3%E3%82%B4%E3%83%A1%E3%83%AA%E4%B9%97%E7%AE%97
// https://www.researchgate.net/publication/3044233_The_Montgomery_modular_inverse_-_Revisited


use std::ops::{Add, Sub, Mul, Div, Rem, Neg, Shl, Shr, AddAssign, SubAssign, MulAssign, DivAssign, RemAssign, ShlAssign, ShrAssign};
use num_traits::{Zero, One, Inv};
use std::collections::BTreeMap;


pub trait Square {
    type Output;
    
    fn square(self) -> Self::Output;
}


const fn inv_mod_power_of_two(a: u64) -> u64 {
    let mut x: u64 = 1;
    let mut i = 0;
    
    while i < u64::BITS.ilog2() {
        x = x.wrapping_mul(2_u64.wrapping_sub(a.wrapping_mul(x)));
        i += 1;
    }
    
    x
}


#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct MontgomeryForm<const N: u64>(u64);

impl<const N: u64> MontgomeryForm<N> {
    const R2: u64 = {
        let r = N.wrapping_neg() as u128;
        ((r * r) % (N as u128)) as u64
    };
    
    const R3: u64 = {
        let r = N.wrapping_neg() as u128;
        let r2 = MontgomeryForm::<N>::R2 as u128;
        ((r2 * r) % (N as u128)) as u64
    };
    
    const N_: u64 = inv_mod_power_of_two(N).wrapping_neg();
    
    pub const ZERO: MontgomeryForm<N> = MontgomeryForm(0);
    
    pub const ONE: MontgomeryForm<N> = {
        let x = MontgomeryForm::<N>::reduce_small(MontgomeryForm::<N>::R2);
        MontgomeryForm(x)
    };
    
    const fn reduce(t: u128) -> u64 {
        debug_assert!(t < (N as u128) << u64::BITS);
        let x = (t as u64).wrapping_mul(MontgomeryForm::<N>::N_);
        let x = (x as u128) * (N as u128);
        let (x, f) = x.overflowing_add(t);
        let x = (x >> u64::BITS) as u64;
        if f || x >= N { x.wrapping_sub(N) } else { x }
    }
    
    const fn reduce_small(t: u64) -> u64 {
        debug_assert!(t < N);
        let x = t.wrapping_mul(MontgomeryForm::<N>::N_);
        let x = (x as u128) * (N as u128);
        let x = x.wrapping_add(t as u128);
        let x = (x >> u64::BITS) as u64;
        x
    }
    
    const fn almost_montgomery_inverse(a: u64) -> (u64, u32) {
        if a == 0 {
            panic!();
        }
        
        let (mut u, mut v) = (N, a);
        let (mut r, mut s): (u64, u64) = (0, 1);
        let mut k = 0;
        loop {
            let n = v.trailing_zeros();
            (v, r, k) = (v >> n, r << n, k + n);
            
            if u > v {
                (u, r, s, k) = ((u - v) / 2, r + s, s * 2, k + 1);
            }
            else if u < v {
                (v, s, r, k) = ((v - u) / 2, s + r, r * 2, k + 1);
            }
            else {
                break;
            }
            
            let n = u.trailing_zeros();
            (u, s, k) = (u >> n, s << n, k + n);
        }
        k += 1;
        let (r, f) = r.overflowing_mul(2);
        let r = if f || r >= N { r.wrapping_sub(N) } else { r };
        let r = N - r;
        
        (r, k)
    }
    
    const fn montgomery_inverse(a: u64) -> u64 {
        let (r, k) = MontgomeryForm::<N>::almost_montgomery_inverse(a);
        
        let m = u64::BITS;
        let (r, e) = if k <= m {
            let r2 = MontgomeryForm::<N>::R2 as u128;
            (MontgomeryForm::<N>::reduce((r as u128) * r2), m - k)
        }
        else {
            (r, 2 * m - k)
        };
        let x = MontgomeryForm::<N>::reduce((r as u128) << e);
        
        x
    }
}

impl<const N: u64> From<u64> for MontgomeryForm<N> {
    fn from(a: u64) -> MontgomeryForm<N> {
        let x = (a as u128) * (MontgomeryForm::<N>::R2 as u128);
        let x = MontgomeryForm::<N>::reduce(x);
        MontgomeryForm(x)
    }
}

impl<const N: u64> From<MontgomeryForm<N>> for u64 {
    fn from(a: MontgomeryForm<N>) -> u64 {
        MontgomeryForm::<N>::reduce_small(a.0)
    }
}

impl<const N: u64> Add<MontgomeryForm<N>> for MontgomeryForm<N> {
    type Output = MontgomeryForm<N>;
    
    fn add(self, other: MontgomeryForm<N>) -> MontgomeryForm<N> {
        let (a, b) = (self.0, other.0);
        let (x, f) = a.overflowing_add(b);
        let x = if f || x >= N { x.wrapping_sub(N) } else { x };
        MontgomeryForm(x)
    }
}

impl<const N: u64> Sub<MontgomeryForm<N>> for MontgomeryForm<N> {
    type Output = MontgomeryForm<N>;
    
    fn sub(self, other: MontgomeryForm<N>) -> MontgomeryForm<N> {
        let (a, b) = (self.0, other.0);
        let (x, f) = a.overflowing_sub(b);
        let x = if f { x.wrapping_add(N) } else { x };
        MontgomeryForm(x)
    }
}

impl<const N: u64> Mul<MontgomeryForm<N>> for MontgomeryForm<N> {
    type Output = MontgomeryForm<N>;
    
    fn mul(self, other: MontgomeryForm<N>) -> MontgomeryForm<N> {
        let (a, b) = (self.0, other.0);
        let x = (a as u128) * (b as u128);
        let x = MontgomeryForm::<N>::reduce(x);
        MontgomeryForm(x)
    }
}

impl<const N: u64> Mul<usize> for MontgomeryForm<N> {
    type Output = MontgomeryForm<N>;
    
    fn mul(self, other: usize) -> MontgomeryForm<N> {
        assert!(usize::BITS <= u64::BITS);
        let b = MontgomeryForm::<N>::from(other as u64);
        self * b
    }
}

impl<const N: u64> Div<MontgomeryForm<N>> for MontgomeryForm<N> {
    type Output = MontgomeryForm<N>;
    
    fn div(self, other: MontgomeryForm<N>) -> MontgomeryForm<N> {
        self * other.inv()
    }
}

impl<const N: u64> Div<usize> for MontgomeryForm<N> {
    type Output = MontgomeryForm<N>;
    
    fn div(self, other: usize) -> MontgomeryForm<N> {
        assert!(usize::BITS <= u64::BITS);
        let x = MontgomeryForm::<N>::montgomery_inverse(other as u64);
        let x = MontgomeryForm(x);
        self * x
    }
}

impl<const N: u64> Neg for MontgomeryForm<N> {
    type Output = MontgomeryForm<N>;
    
    fn neg(self) -> MontgomeryForm<N> {
        let a = self.0;
        let (x, f) = a.overflowing_neg();
        let x = if f { x.wrapping_add(N) } else { x };
        MontgomeryForm(x)
    }
}

impl<const N: u64> AddAssign<MontgomeryForm<N>> for MontgomeryForm<N> {
    fn add_assign(&mut self, other: MontgomeryForm<N>) {
        *self = *self + other;
    }
}

impl<const N: u64> SubAssign<MontgomeryForm<N>> for MontgomeryForm<N> {
    fn sub_assign(&mut self, other: MontgomeryForm<N>) {
        *self = *self - other;
    }
}

impl<const N: u64> MulAssign<MontgomeryForm<N>> for MontgomeryForm<N> {
    fn mul_assign(&mut self, other: MontgomeryForm<N>) {
        *self = *self * other;
    }
}

impl<const N: u64> DivAssign<MontgomeryForm<N>> for MontgomeryForm<N> {
    fn div_assign(&mut self, other: MontgomeryForm<N>) {
        *self = *self / other;
    }
}

impl<const N: u64> Zero for MontgomeryForm<N> {
    fn zero() -> MontgomeryForm<N> {
        MontgomeryForm::<N>::ZERO
    }
    
    fn is_zero(&self) -> bool {
        *self == MontgomeryForm::<N>::ZERO
    }
}

impl<const N: u64> One for MontgomeryForm<N> {
    fn one() -> MontgomeryForm<N> {
        MontgomeryForm::<N>::ONE
    }
    
    fn is_one(&self) -> bool {
        *self == MontgomeryForm::<N>::ONE
    }
}

impl<const N: u64> Inv for MontgomeryForm<N> {
    type Output = MontgomeryForm<N>;
    
    fn inv(self) -> MontgomeryForm<N> {
        let (r, k) = MontgomeryForm::<N>::almost_montgomery_inverse(self.0);
        
        let m = u64::BITS;
        let (c, e) = if k <= m {
            (MontgomeryForm::<N>::R3 as u128, m - k)
        }
        else {
            (MontgomeryForm::<N>::R2 as u128, 2 * m - k)
        };
        let r = MontgomeryForm::<N>::reduce((r as u128) * c);
        let x = MontgomeryForm::<N>::reduce((r as u128) << e);
        
        MontgomeryForm(x)
    }
}


const SSA_THRESHOLD: usize = 100;


trait Coefficient
where
    Self: Copy + PartialEq + Zero + One + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Neg<Output = Self> + Inv<Output = Self> + Mul<usize, Output = Self> + Div<usize, Output = Self> + AddAssign + SubAssign + MulAssign,
{
}

impl<T> Coefficient for T
where
    T: Copy + PartialEq + Zero + One + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Neg<Output = T> + Inv<Output = T> + Mul<usize, Output = T> + Div<usize, Output = T> + AddAssign + SubAssign + MulAssign,
{
}


#[derive(Clone, PartialEq, Eq, Debug)]
struct Polynomial<T>(Vec<T>);

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct PolynomialRef<'a, T>(&'a [T]);

impl<T> Polynomial<T>
where
    T: Coefficient,
{
    pub fn as_borrowed_polynomial(&self) -> PolynomialRef<'_, T> {
        PolynomialRef(self.0.as_slice())
    }
    
    pub fn len(&self) -> usize {
        self.as_borrowed_polynomial().len()
    }
    
    pub fn degree(&self) -> Option<usize> {
        self.as_borrowed_polynomial().degree()
    }
    
    pub fn coefficient(&self, e: usize) -> T {
        self.as_borrowed_polynomial().coefficient(e)
    }
    
    pub fn x() -> Polynomial<T> {
        Polynomial(vec![T::zero(), T::one()])
    }
    
    fn div_power_of_x(self, e: usize) -> Polynomial<T> {
        if e >= self.len() {
            return Polynomial::zero();
        }
        
        let mut v = self.0;
        v.copy_within(e.., 0);
        v.truncate(v.len() - e);
        Polynomial(v)
    }
    
    fn mod_power_of_x(self, e: usize) -> Polynomial<T> {
        if e >= self.len() {
            return self;
        }
        
        let mut v = self.0;
        v.truncate(e);
        Polynomial::from(v)
    }
    
    pub fn monic(self) -> Polynomial<T> {
        let mut a = self;
        
        let Some(&l) = a.0.last() else {
            return Polynomial::zero();
        };
        let c = l.inv();
        for x in a.0.iter_mut() {
            *x *= c;
        }
        
        a
    }
    
    fn inv_mod_power_of_x(self, e: usize) -> Polynomial<T> {
        self.as_borrowed_polynomial().inv_mod_power_of_x(e)
    }
    
    fn rev_mod_power_of_x(self, k: usize, n: usize) -> Polynomial<T> {
        let Some(d) = self.degree() else {
            return Polynomial::zero();
        };
        if k < d {
            panic!();
        }
        let e = k - d;
        if e >= n {
            return Polynomial::zero();
        }
        
        let l1 = std::cmp::min(k + 1, n);
        if self.0.capacity() <= l1 {
            return self.as_borrowed_polynomial().rev_mod_power_of_x(k, n);
        }
        
        let mut v = self.0;
        let l2 = std::cmp::max(l1, d + 1);
        v.resize(l2, T::zero());
        v[k + 1 - l2..].reverse();
        v[..e].fill(T::zero());
        v.truncate(l1);
        
        Polynomial::from(v)
    }
    
    pub fn div_rem(&self, other: &Polynomial<T>) -> (Polynomial<T>, Polynomial<T>) {
        self.as_borrowed_polynomial().div_rem(other.as_borrowed_polynomial())
    }
    
    fn half_gcd(&self, other: &Polynomial<T>) -> (Polynomial<T>, Polynomial<T>, Polynomial<T>, Polynomial<T>) {
        self.as_borrowed_polynomial().half_gcd(other.as_borrowed_polynomial())
    }
    
    pub fn gcd(&self, other: &Polynomial<T>) -> Polynomial<T> {
        self.as_borrowed_polynomial().gcd(other.as_borrowed_polynomial())
    }
}

impl<'a, T> PolynomialRef<'a, T>
where
    T: Coefficient,
{
    pub fn to_owned_polynomial(self) -> Polynomial<T> {
        Polynomial(self.0.to_vec())
    }
    
    pub fn len(self) -> usize {
        self.0.len()
    }
    
    pub fn degree(self) -> Option<usize> {
        self.len().checked_sub(1)
    }
    
    pub fn coefficient(self, e: usize) -> T {
        self.0.get(e).copied().unwrap_or(T::zero())
    }
    
    fn div_power_of_x(self, e: usize) -> PolynomialRef<'a, T> {
        if e >= self.len() {
            return PolynomialRef(&[]);
        }
        
        PolynomialRef(&self.0[e..])
    }
    
    fn mod_power_of_x(self, e: usize) -> PolynomialRef<'a, T> {
        if e >= self.len() {
            return self;
        }
        
        PolynomialRef::from(&self.0[..e])
    }
    
    fn inv_mod_power_of_x(self, e: usize) -> Polynomial<T> {
        let f = self;
        let Some(t) = f.0.first() else {
            panic!();
        };
        if e == 0 {
            return Polynomial::zero();
        }
        
        let mut g = Polynomial(vec![t.inv()]);
        if e <= 1 {
            return g;
        }
        
        for i in (0..=(e - 1).ilog2()).rev() {
            let n = e.div_ceil(1 << i);
            let f2 = f.mod_power_of_x(n);
            let r = (&g).square().mod_power_of_x(n);
            let r = (f2 * r).mod_power_of_x(n);
            
            let g2 = g * 2 - r;
            g = g2;
        }
        
        g
    }
    
    fn rev_mod_power_of_x(self, k: usize, n: usize) -> Polynomial<T> {
        let Some(d) = self.degree() else {
            return Polynomial::zero();
        };
        if k < d {
            panic!();
        }
        let e = k - d;
        if e >= n {
            return Polynomial::zero();
        }
        
        let mut v = Vec::with_capacity(std::cmp::min(k + 1, n));
        v.resize(e, T::zero());
        v.extend(self.0.iter().copied().rev().take(n - e));
        return Polynomial::from(v)
    }
    
    pub fn div_rem(self, other: PolynomialRef<'_, T>) -> (Polynomial<T>, Polynomial<T>) {
        let (a, b) = (self, other);
        if a.len() < b.len() {
            return (Polynomial::zero(), a.to_owned_polynomial());
        }
        let Some(m) = b.degree() else {
            panic!();
        };
        
        let q = a / b;
        let q2 = q.as_borrowed_polynomial().mod_power_of_x(m);
        let r = a.mod_power_of_x(m) - (b * q2).mod_power_of_x(m);
        
        (q, r)
    }
    
    fn half_gcd(self, other: PolynomialRef<'_, T>) -> (Polynomial<T>, Polynomial<T>, Polynomial<T>, Polynomial<T>) {
        let (a, b) = (self, other);
        
        let Some(d) = a.degree() else {
            return (Polynomial::one(), Polynomial::zero(), Polynomial::zero(), Polynomial::one());
        };
        let m = d.div_ceil(2);
        if !b.degree().is_some_and(|x| x >= m) {
            return (Polynomial::one(), Polynomial::zero(), Polynomial::zero(), Polynomial::one());
        }
        
        let a0 = a.div_power_of_x(m);
        let b0 = b.div_power_of_x(m);
        let (r11, r12, r21, r22) = a0.half_gcd(b0);
        let a_ = &r11 * a + &r12 * b;
        let b_ = &r21 * a + &r22 * b;
        let Some(l) = b_.degree().filter(|&x| x >= m) else {
            return (r11, r12, r21, r22);
        };
        
        let (q, d) = a_.div_rem(&b_);
        let c = b_;
        let t21 = r11 - &q * &r21;
        let t22 = r12 - &q * &r22;
        let (t11, t12) = (r21, r22);
        if !d.degree().is_some_and(|x| x >= m) {
            return (t11, t12, t21, t22);
        }
        
        let k = 2 * m - l;
        let c0 = c.div_power_of_x(k);
        let d0 = d.div_power_of_x(k);
        let (s11, s12, s21, s22) = c0.half_gcd(&d0);
        let m11 = &s11 * &t11 + &s12 * &t21;
        let m12 = &s11 * &t12 + &s12 * &t22;
        let m21 = &s21 * &t11 + &s22 * &t21;
        let m22 = &s21 * &t12 + &s22 * &t22;
        (m11, m12, m21, m22)
    }
    
    pub fn gcd(self, other: PolynomialRef<'_, T>) -> Polynomial<T> {
        let (a, b) = (self, other);
        let (a, b) = if a.len() >= b.len() { (a, b) } else { (b, a) };
        if b.degree().is_none() {
            return a.to_owned_polynomial();
        }
        
        let (m11, m12, m21, m22) = a.half_gcd(b);
        let a_ = m11 * a + m12 * b;
        let b_ = m21 * a + m22 * b;
        if b_.degree().is_none() {
            return a_;
        }
        
        let d = a_ % &b_;
        let c = b_;
        c.gcd(&d)
    }
}

impl<T> From<Vec<T>> for Polynomial<T>
where
    T: Coefficient,
{
    fn from(x: Vec<T>) -> Polynomial<T> {
        let mut x = x;
        let n = x.iter().rev().skip_while(|c| c.is_zero()).count();
        x.truncate(n);
        Polynomial(x)
    }
}

impl<T> From<Polynomial<T>> for Vec<T>
where
    T: Coefficient,
{
    fn from(x: Polynomial<T>) -> Vec<T> {
        x.0
    }
}

impl<'a, T> From<&'a [T]> for PolynomialRef<'a, T>
where
    T: Coefficient,
{
    fn from(x: &[T]) -> PolynomialRef<'_, T> {
        let n = x.iter().rev().skip_while(|x| x.is_zero()).count();
        PolynomialRef(&x[..n])
    }
}

impl<'a, T> From<PolynomialRef<'a, T>> for &'a [T]
where
    T: Coefficient,
{
    fn from(x: PolynomialRef<'_, T>) -> &[T] {
        x.0
    }
}

impl<T> Add<PolynomialRef<'_, T>> for PolynomialRef<'_, T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn add(self, other: PolynomialRef<'_, T>) -> Polynomial<T> {
        let (a, b) = (self.0, other.0);
        let (a, b) = if a.len() >= b.len() { (a, b) } else { (b, a) };
        
        let mut v = Vec::with_capacity(a.len());
        let (a1, a2) = a.split_at(b.len());
        for (&x, &y) in std::iter::zip(a1.iter(), b.iter()) {
            v.push(x + y);
        }
        v.extend_from_slice(a2);
        
        Polynomial::from(v)
    }
}

impl<T> Sub<PolynomialRef<'_, T>> for PolynomialRef<'_, T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn sub(self, other: PolynomialRef<'_, T>) -> Polynomial<T> {
        let (a, b) = (self.0, other.0);
        
        if a.len() >= b.len() {
            let mut v = Vec::with_capacity(a.len());
            let (a1, a2) = a.split_at(b.len());
            for (&x, &y) in std::iter::zip(a1.iter(), b.iter()) {
                v.push(x - y);
            }
            v.extend_from_slice(a2);
            Polynomial::from(v)
        }
        else {
            let mut v = Vec::with_capacity(b.len());
            let (b1, b2) = b.split_at(a.len());
            for (&x, &y) in std::iter::zip(a.iter(), b1.iter()) {
                v.push(x - y);
            }
            for &y in b2.iter() {
                v.push(-y);
            }
            Polynomial::from(v)
        }
    }
}

impl<T> Mul<PolynomialRef<'_, T>> for PolynomialRef<'_, T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn mul(self, other: PolynomialRef<'_, T>) -> Polynomial<T> {
        let (a, b) = (self.0, other.0);
        let (a, b) = if a.len() <= b.len() { (a, b) } else { (b, a) };
        if a.len() == 0 {
            return Polynomial::zero();
        }
        
        if a.len() < SSA_THRESHOLD {
            let r_len = a.len() + b.len() - 1;
            let mut r = vec![T::zero(); r_len];
            for (i, &x) in a.iter().enumerate() {
                for (&y, z) in std::iter::zip(b.iter(), r[i..i + b.len()].iter_mut()) {
                    *z += x * y;
                }
            }
            return Polynomial(r);
        }
        
        let r = ssa_mul(a, b);
        Polynomial(r)
    }
}

impl<T> Div<PolynomialRef<'_, T>> for PolynomialRef<'_, T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn div(self, other: PolynomialRef<'_, T>) -> Polynomial<T> {
        let (a, b) = (self, other);
        if a.len() < b.len() {
            return Polynomial::zero();
        }
        let (Some(n), Some(m)) = (a.degree(), b.degree()) else {
            panic!();
        };
        
        let h = n - m + 1;
        let f = b.rev_mod_power_of_x(m, h);
        let g = f.inv_mod_power_of_x(h);
        let a2 = a.rev_mod_power_of_x(n, h);
        let q = (a2 * g).mod_power_of_x(h);
        let q = q.rev_mod_power_of_x(n - m, h);
        
        q
    }
}

impl<T> Rem<PolynomialRef<'_, T>> for PolynomialRef<'_, T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn rem(self, other: PolynomialRef<'_, T>) -> Polynomial<T> {
        let (_, r) = self.div_rem(other);
        r
    }
}

macro_rules! impl_polynomial_binop {
    (
        $Trait:ident, $method:ident, $TraitAssign:ident, $method_assign:ident
    ) => {
        impl<T> $Trait<Polynomial<T>> for Polynomial<T>
        where
            T: Coefficient,
        {
            type Output = Polynomial<T>;
            
            fn $method(self, other: Polynomial<T>) -> Polynomial<T> {
                self.as_borrowed_polynomial().$method(other.as_borrowed_polynomial())
            }
        }
        
        impl<T> $Trait<&Polynomial<T>> for Polynomial<T>
        where
            T: Coefficient,
        {
            type Output = Polynomial<T>;
            
            fn $method(self, other: &Polynomial<T>) -> Polynomial<T> {
                self.as_borrowed_polynomial().$method(other.as_borrowed_polynomial())
            }
        }
        
        impl<T> $Trait<Polynomial<T>> for &Polynomial<T>
        where
            T: Coefficient,
        {
            type Output = Polynomial<T>;
            
            fn $method(self, other: Polynomial<T>) -> Polynomial<T> {
                self.as_borrowed_polynomial().$method(other.as_borrowed_polynomial())
            }
        }
        
        impl<T> $Trait<&Polynomial<T>> for &Polynomial<T>
        where
            T: Coefficient,
        {
            type Output = Polynomial<T>;
            
            fn $method(self, other: &Polynomial<T>) -> Polynomial<T> {
                self.as_borrowed_polynomial().$method(other.as_borrowed_polynomial())
            }
        }
        
        impl<T> $Trait<PolynomialRef<'_, T>> for Polynomial<T>
        where
            T: Coefficient,
        {
            type Output = Polynomial<T>;
            
            fn $method(self, other: PolynomialRef<'_, T>) -> Polynomial<T> {
                self.as_borrowed_polynomial().$method(other)
            }
        }
        
        impl<T> $Trait<PolynomialRef<'_, T>> for &Polynomial<T>
        where
            T: Coefficient,
        {
            type Output = Polynomial<T>;
            
            fn $method(self, other: PolynomialRef<'_, T>) -> Polynomial<T> {
                self.as_borrowed_polynomial().$method(other)
            }
        }
        
        impl<T> $Trait<Polynomial<T>> for PolynomialRef<'_, T>
        where
            T: Coefficient,
        {
            type Output = Polynomial<T>;
            
            fn $method(self, other: Polynomial<T>) -> Polynomial<T> {
                self.$method(other.as_borrowed_polynomial())
            }
        }
        
        impl<T> $Trait<&Polynomial<T>> for PolynomialRef<'_, T>
        where
            T: Coefficient,
        {
            type Output = Polynomial<T>;
            
            fn $method(self, other: &Polynomial<T>) -> Polynomial<T> {
                self.$method(other.as_borrowed_polynomial())
            }
        }
        
        impl<T> $TraitAssign<Polynomial<T>> for Polynomial<T>
        where
            T: Coefficient,
        {
            fn $method_assign(&mut self, other: Polynomial<T>) {
                *self = self.as_borrowed_polynomial().$method(other.as_borrowed_polynomial());
            }
        }
        
        impl<T> $TraitAssign<&Polynomial<T>> for Polynomial<T>
        where
            T: Coefficient,
        {
            fn $method_assign(&mut self, other: &Polynomial<T>) {
                *self = self.as_borrowed_polynomial().$method(other.as_borrowed_polynomial());
            }
        }
        
        impl<T> $TraitAssign<PolynomialRef<'_, T>> for Polynomial<T>
        where
            T: Coefficient,
        {
            fn $method_assign(&mut self, other: PolynomialRef<'_, T>) {
                *self = self.as_borrowed_polynomial().$method(other);
            }
        }
    };
}

impl_polynomial_binop! {
    Add, add, AddAssign, add_assign
}

impl_polynomial_binop! {
    Sub, sub, SubAssign, sub_assign
}

impl_polynomial_binop! {
    Mul, mul, MulAssign, mul_assign
}

impl_polynomial_binop! {
    Div, div, DivAssign, div_assign
}

impl_polynomial_binop! {
    Rem, rem, RemAssign, rem_assign
}

impl<T> Mul<usize> for Polynomial<T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn mul(self, other: usize) -> Polynomial<T> {
        let mut a = self;
        let c = T::one() * other;
        for x in a.0.iter_mut() {
            *x *= c;
        }
        a
    }
}

impl<T> Mul<usize> for &Polynomial<T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn mul(self, other: usize) -> Polynomial<T> {
        self.as_borrowed_polynomial().mul(other)
    }
}

impl<T> Mul<usize> for PolynomialRef<'_, T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn mul(self, other: usize) -> Polynomial<T> {
        self.to_owned_polynomial().mul(other)
    }
}

impl<T> MulAssign<usize> for Polynomial<T>
where
    T: Coefficient,
{
    fn mul_assign(&mut self, other: usize) {
        let v = std::mem::take(&mut self.0);
        let (a, b) = (Polynomial(v), other);
        *self = a.mul(b);
    }
}

impl<T> Shl<usize> for Polynomial<T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn shl(self, other: usize) -> Polynomial<T> {
        let (a, e) = (self, other);
        if a.0.capacity() - a.0.len() < e {
            return a.as_borrowed_polynomial().shl(e);
        }
        
        let mut v = a.0;
        let n = v.len();
        if e >= n {
            v.resize(e, T::zero());
            v.extend_from_within(..n);
            v[..n].fill(T::zero());
        }
        else {
            v.extend_from_within(n - e..n);
            v.copy_within(..n - e, e);
            v[..e].fill(T::zero());
        }
        
        Polynomial(v)
    }
}

impl<T> Shl<usize> for &Polynomial<T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn shl(self, other: usize) -> Polynomial<T> {
        self.as_borrowed_polynomial().shl(other)
    }
}

impl<T> Shl<usize> for PolynomialRef<'_, T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn shl(self, other: usize) -> Polynomial<T> {
        let (a, e) = (self, other);
        
        let mut v = Vec::with_capacity(a.len() + e);
        v.resize(e, T::zero());
        v.extend_from_slice(a.0);
        Polynomial(v)
    }
}

impl<T> ShlAssign<usize> for Polynomial<T>
where
    T: Coefficient,
{
    fn shl_assign(&mut self, other: usize) {
        let v = std::mem::take(&mut self.0);
        let (a, b) = (Polynomial(v), other);
        *self = a.shl(b);
    }
}

impl<T> Shr<usize> for Polynomial<T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn shr(self, other: usize) -> Polynomial<T> {
        self.div_power_of_x(other)
    }
}

impl<T> Shr<usize> for &Polynomial<T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn shr(self, other: usize) -> Polynomial<T> {
        self.as_borrowed_polynomial().shr(other)
    }
}

impl<T> Shr<usize> for PolynomialRef<'_, T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn shr(self, other: usize) -> Polynomial<T> {
        self.div_power_of_x(other).to_owned_polynomial()
    }
}

impl<T> ShrAssign<usize> for Polynomial<T>
where
    T: Coefficient,
{
    fn shr_assign(&mut self, other: usize) {
        let v = std::mem::take(&mut self.0);
        let (a, b) = (Polynomial(v), other);
        *self = a.shr(b);
    }
}

impl<T> Zero for Polynomial<T>
where
    T: Coefficient,
{
    fn zero() -> Polynomial<T> {
        Polynomial(vec![])
    }
    
    fn is_zero(&self) -> bool {
        self.len() == 0
    }
}

impl<T> One for Polynomial<T>
where
    T: Coefficient,
{
    fn one() -> Polynomial<T> {
        Polynomial(vec![T::one()])
    }
    
    fn is_one(&self) -> bool {
        matches!(*self.0, [x] if x.is_one())
    }
}

impl<T> Square for Polynomial<T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn square(self) -> Polynomial<T> {
        self.as_borrowed_polynomial().square()
    }
}

impl<T> Square for &Polynomial<T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn square(self) -> Polynomial<T> {
        self.as_borrowed_polynomial().square()
    }
}

impl<T> Square for PolynomialRef<'_, T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn square(self) -> Polynomial<T> {
        let a = self.0;
        if a.len() == 0 {
            return Polynomial::zero();
        }
        
        if a.len() < SSA_THRESHOLD {
            let r_len = 2 * a.len() - 1;
            let mut r = vec![T::zero(); r_len];
            let c = T::one() * 2;
            for (i, &x) in a.iter().enumerate() {
                let a2 = &a[..i];
                for (&y, z) in std::iter::zip(a2.iter(), r[i..i + a2.len()].iter_mut()) {
                    *z += c * x * y;
                }
                let z = &mut r[2 * i];
                *z += x * x;
            }
            return Polynomial(r);
        }
        
        let r = ssa_square(a);
        Polynomial(r)
    }
}


fn ntt<T>(buf: &mut [T], w: usize, n: usize, k: usize)
where
    T: Coefficient,
{
    debug_assert!(w.is_power_of_two() && n.is_power_of_two() && k.is_power_of_two());
    debug_assert_eq!(w * k, 2 * n);
    debug_assert_eq!(buf.len(), k * n);
    if k <= 1 {
        return;
    }
    
    let mut temp_buf = vec![T::zero(); n];
    ntt_recursion(buf, w, n, k, &mut temp_buf);
}

fn ntt_recursion<T>(buf: &mut [T], w: usize, n: usize, k: usize, temp_buf: &mut [T])
where
    T: Coefficient,
{
    let k1 = k / 2;
    let (buf_1, buf_2) = buf.split_at_mut(k1 * n);
    let r = &mut *temp_buf;
    for (i, (a, b)) in std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks_mut(n)).enumerate() {
        for ((&x, &y), z) in std::iter::zip(std::iter::zip(a.iter(), b.iter()), r.iter_mut()) {
            *z = x - y;
        }
        for (x, &y) in std::iter::zip(a.iter_mut(), b.iter()) {
            *x += y;
        }
        let (b1, b2) = b.split_at_mut(w * i);
        let (r1, r2) = r.split_at(n - w * i);
        b2.copy_from_slice(r1);
        for (x, &y) in std::iter::zip(b1.iter_mut(), r2.iter()) {
            *x = -y;
        }
    }
    if k1 > 1 {
        ntt_recursion(buf_1, 2 * w, n, k1, temp_buf);
        ntt_recursion(buf_2, 2 * w, n, k1, temp_buf);
    }
}


fn intt<T>(buf: &mut [T], w: usize, n: usize, k: usize)
where
    T: Coefficient,
{
    debug_assert!(w.is_power_of_two() && n.is_power_of_two() && k.is_power_of_two());
    debug_assert_eq!(w * k, 2 * n);
    debug_assert_eq!(buf.len(), k * n);
    if k <= 1 {
        return;
    }
    
    let mut temp_buf = vec![T::zero(); n];
    intt_recursion(buf, w, n, k, &mut temp_buf);
    
    let c = T::one() / k;
    for x in buf.iter_mut() {
        *x *= c;
    }
}

fn intt_recursion<T>(buf: &mut [T], w: usize, n: usize, k: usize, temp_buf: &mut [T])
where
    T: Coefficient,
{
    let k1 = k / 2;
    let (buf_1, buf_2) = buf.split_at_mut(k1 * n);
    if k1 > 1 {
        intt_recursion(buf_1, 2 * w, n, k1, temp_buf);
        intt_recursion(buf_2, 2 * w, n, k1, temp_buf);
    }
    let r = &mut *temp_buf;
    for (i, (a, b)) in std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks_mut(n)).enumerate() {
        let (b1, b2) = b.split_at(w * i);
        let (r1, r2) = r.split_at_mut(n - w * i);
        r1.copy_from_slice(b2);
        for (x, &y) in std::iter::zip(r2.iter_mut(), b1.iter()) {
            *x = -y;
        }
        for ((&x, &y), z) in std::iter::zip(std::iter::zip(a.iter(), r.iter()), b.iter_mut()) {
            *z = x - y;
        }
        for (x, &y) in std::iter::zip(a.iter_mut(), r.iter()) {
            *x += y;
        }
    }
}


fn truncated_ntt<T>(buf: &mut [T], w: usize, n: usize, k: usize, k2: usize, k3: usize)
where
    T: Coefficient,
{
    debug_assert!(w.is_power_of_two() && n.is_power_of_two() && k.is_power_of_two());
    debug_assert_eq!(w * k, 2 * n);
    debug_assert_eq!(buf.len(), k * n);
    debug_assert!(k2 <= k && k3 <= k2);
    if k <= 1 || k2 == 0 {
        return;
    }
    
    let mut w = w;
    let mut k = k;
    let mut k1 = k / 2;
    while k2 <= k1 {
        k /= 2;
        k1 /= 2;
        w *= 2;
    }
    let buf = &mut buf[..k * n];
    
    let mut temp_buf = vec![T::zero(); n];
    let (buf_1, buf_2) = buf.split_at_mut(k1 * n);
    if k3 <= k1 {
        let (buf_11, buf_12) = buf_1.split_at_mut(k3 * n);
        let (buf_21, buf_22) = buf_2.split_at_mut(k3 * n);
        for (i, (a, b)) in std::iter::zip(buf_11.chunks(n), buf_21.chunks_mut(n)).enumerate() {
            let (b1, b2) = b.split_at_mut(w * i);
            let (a1, a2) = a.split_at(n - w * i);
            b2.copy_from_slice(a1);
            for (x, &y) in std::iter::zip(b1.iter_mut(), a2.iter()) {
                *x = -y;
            }
        }
        buf_12.fill(T::zero());
        buf_22.fill(T::zero());
        truncated_ntt_recursion_2(buf_1, 2 * w, n, k1, k2, k3, &mut temp_buf);
        truncated_ntt_recursion_2(buf_2, 2 * w, n, k1, k2 - k1, k3, &mut temp_buf);
    }
    else {
        let r = &mut temp_buf;
        let mut it = std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks_mut(n)).enumerate();
        for (i, (a, b)) in it.by_ref().take(k3) {
            for ((&x, &y), z) in std::iter::zip(std::iter::zip(a.iter(), b.iter()), r.iter_mut()) {
                *z = x - y;
            }
            for (x, &y) in std::iter::zip(a.iter_mut(), b.iter()) {
                *x += y;
            }
            let (b1, b2) = b.split_at_mut(w * i);
            let (r1, r2) = r.split_at(n - w * i);
            b2.copy_from_slice(r1);
            for (x, &y) in std::iter::zip(b1.iter_mut(), r2.iter()) {
                *x = -y;
            }
        }
        for (i, (a, b)) in it {
            let (b1, b2) = b.split_at_mut(w * i);
            let (a1, a2) = a.split_at(n - w * i);
            b2.copy_from_slice(a1);
            for (x, &y) in std::iter::zip(b1.iter_mut(), a2.iter()) {
                *x = -y;
            }
        }
        truncated_ntt_recursion(buf_1, 2 * w, n, k1, k2, &mut temp_buf);
        truncated_ntt_recursion(buf_2, 2 * w, n, k1, k2 - k1, &mut temp_buf);
    }
}

fn truncated_ntt_recursion<T>(buf: &mut [T], w: usize, n: usize, k: usize, k2: usize, temp_buf: &mut [T])
where
    T: Coefficient,
{
    if k <= 1 {
        return;
    }
    let k1 = k / 2;
    let (buf_1, buf_2) = buf.split_at_mut(k1 * n);
    if k2 <= k1 {
        for (a, b) in std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks(n)) {
            for (x, &y) in std::iter::zip(a.iter_mut(), b.iter()) {
                *x += y;
            }
        }
        truncated_ntt_recursion(buf_1, 2 * w, n, k1, k2, temp_buf);
    }
    else {
        let r = &mut *temp_buf;
        for (i, (a, b)) in std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks_mut(n)).enumerate() {
            for ((&x, &y), z) in std::iter::zip(std::iter::zip(a.iter(), b.iter()), r.iter_mut()) {
                *z = x - y;
            }
            for (x, &y) in std::iter::zip(a.iter_mut(), b.iter()) {
                *x += y;
            }
            let (b1, b2) = b.split_at_mut(w * i);
            let (r1, r2) = r.split_at(n - w * i);
            b2.copy_from_slice(r1);
            for (x, &y) in std::iter::zip(b1.iter_mut(), r2.iter()) {
                *x = -y;
            }
        }
        truncated_ntt_recursion(buf_1, 2 * w, n, k1, k2, temp_buf);
        truncated_ntt_recursion(buf_2, 2 * w, n, k1, k2 - k1, temp_buf);
    }
}

fn truncated_ntt_recursion_2<T>(buf: &mut [T], w: usize, n: usize, k: usize, k2: usize, k3: usize, temp_buf: &mut [T])
where
    T: Coefficient,
{
    if k <= 1 {
        return;
    }
    let k1 = k / 2;
    let (buf_1, buf_2) = buf.split_at_mut(k1 * n);
    if k2 <= k1 {
        if k3 <= k1 {
            truncated_ntt_recursion_2(buf_1, 2 * w, n, k1, k2, k3, temp_buf);
        }
        else {
            for (a, b) in std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks(n)).take(k3) {
                for (x, &y) in std::iter::zip(a.iter_mut(), b.iter()) {
                    *x += y;
                }
            }
            truncated_ntt_recursion(buf_1, 2 * w, n, k1, k2, temp_buf);
        }
    }
    else {
        if k3 <= k1 {
            for (i, (a, b)) in std::iter::zip(buf_1.chunks(n), buf_2.chunks_mut(n)).enumerate().take(k3) {
                let (b1, b2) = b.split_at_mut(w * i);
                let (a1, a2) = a.split_at(n - w * i);
                b2.copy_from_slice(a1);
                for (x, &y) in std::iter::zip(b1.iter_mut(), a2.iter()) {
                    *x = -y;
                }
            }
            truncated_ntt_recursion_2(buf_1, 2 * w, n, k1, k2, k3, temp_buf);
            truncated_ntt_recursion_2(buf_2, 2 * w, n, k1, k2 - k1, k3, temp_buf);
        }
        else {
            let r = &mut *temp_buf;
            let mut it = std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks_mut(n)).enumerate();
            for (i, (a, b)) in it.by_ref().take(k3) {
                for ((&x, &y), z) in std::iter::zip(std::iter::zip(a.iter(), b.iter()), r.iter_mut()) {
                    *z = x - y;
                }
                for (x, &y) in std::iter::zip(a.iter_mut(), b.iter()) {
                    *x += y;
                }
                let (b1, b2) = b.split_at_mut(w * i);
                let (r1, r2) = r.split_at(n - w * i);
                b2.copy_from_slice(r1);
                for (x, &y) in std::iter::zip(b1.iter_mut(), r2.iter()) {
                    *x = -y;
                }
            }
            for (i, (a, b)) in it {
                let (b1, b2) = b.split_at_mut(w * i);
                let (a1, a2) = a.split_at(n - w * i);
                b2.copy_from_slice(a1);
                for (x, &y) in std::iter::zip(b1.iter_mut(), a2.iter()) {
                    *x = -y;
                }
            }
            truncated_ntt_recursion(buf_1, 2 * w, n, k1, k2, temp_buf);
            truncated_ntt_recursion(buf_2, 2 * w, n, k1, k2 - k1, temp_buf);
        }
    }
}


fn truncated_intt<T>(buf: &mut [T], w: usize, n: usize, k: usize, k2: usize)
where
    T: Coefficient,
{
    debug_assert!(w.is_power_of_two() && n.is_power_of_two() && k.is_power_of_two());
    debug_assert_eq!(w * k, 2 * n);
    debug_assert_eq!(buf.len(), k * n);
    debug_assert!(k2 <= k);
    if k <= 1 || k2 == 0 {
        return;
    }
    
    let mut w = w;
    let mut k = k;
    let mut k1 = k / 2;
    while k2 <= k1 {
        k /= 2;
        k1 /= 2;
        w *= 2;
    }
    let buf = &mut buf[..k * n];
    
    let mut temp_buf = vec![T::zero(); n];
    let (buf_1, buf_2) = buf.split_at_mut(k1 * n);
    let c1 = T::one() * 2;
    let c2 = T::one() / 2;
    truncated_intt_recursion(buf_1, 2 * w, n, k1, k2, c1, c2, &mut temp_buf);
    for (i, (a, b)) in std::iter::zip(buf_1.chunks(n), buf_2.chunks_mut(n)).enumerate().skip(k2 - k1) {
        let (b1, b2) = b.split_at_mut(w * i);
        let (a1, a2) = a.split_at(n - w * i);
        b2.copy_from_slice(a1);
        for (x, &y) in std::iter::zip(b1.iter_mut(), a2.iter()) {
            *x = -y;
        }
    }
    truncated_intt_recursion(buf_2, 2 * w, n, k1, k2 - k1, c1, c2, &mut temp_buf);
    let r = &mut temp_buf;
    let mut it = std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks_mut(n)).enumerate();
    for (i, (a, b)) in it.by_ref().take(k2 - k1) {
        let (b1, b2) = b.split_at(w * i);
        let (r1, r2) = r.split_at_mut(n - w * i);
        r1.copy_from_slice(b2);
        for (x, &y) in std::iter::zip(r2.iter_mut(), b1.iter()) {
            *x = -y;
        }
        for ((&x, &y), z) in std::iter::zip(std::iter::zip(a.iter(), r.iter()), b.iter_mut()) {
            *z = c2 * (x - y);
        }
        for (x, &y) in std::iter::zip(a.iter_mut(), r.iter()) {
            *x = c2 * (*x + y);
        }
    }
    for (i, (a, b)) in it {
        let (b1, b2) = b.split_at(w * i);
        let (r1, r2) = r.split_at_mut(n - w * i);
        r1.copy_from_slice(b2);
        for (x, &y) in std::iter::zip(r2.iter_mut(), b1.iter()) {
            *x = -y;
        }
        for (x, &y) in std::iter::zip(a.iter_mut(), r.iter()) {
            *x = c2 * (*x + y);
        }
    }
}

fn truncated_intt_recursion<T>(buf: &mut [T], w: usize, n: usize, k: usize, k2: usize, c1: T, c2: T, temp_buf: &mut [T])
where
    T: Coefficient,
{
    if k <= 1 {
        return;
    }
    let k1 = k / 2;
    let (buf_1, buf_2) = buf.split_at_mut(k1 * n);
    if k2 <= k1 {
        for (a, b) in std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks(n)).skip(k2) {
            for (x, &y) in std::iter::zip(a.iter_mut(), b.iter()) {
                *x += y;
            }
        }
        truncated_intt_recursion(buf_1, 2 * w, n, k1, k2, c1, c2, temp_buf);
        for (a, b) in std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks(n)) {
            for (x, &y) in std::iter::zip(a.iter_mut(), b.iter()) {
                *x -= y;
            }
        }
    }
    else {
        truncated_intt_recursion(buf_1, 2 * w, n, k1, k2, c1, c2, temp_buf);
        let r = &mut *temp_buf;
        for (i, (a, b)) in std::iter::zip(buf_1.chunks(n), buf_2.chunks_mut(n)).enumerate().skip(k2 - k1) {
            for ((&x, &y), z) in std::iter::zip(std::iter::zip(a.iter(), b.iter()), r.iter_mut()) {
                *z = x - c1 * y;
            }
            let (b1, b2) = b.split_at_mut(w * i);
            let (r1, r2) = r.split_at(n - w * i);
            b2.copy_from_slice(r1);
            for (x, &y) in std::iter::zip(b1.iter_mut(), r2.iter()) {
                *x = -y;
            }
        }
        truncated_intt_recursion(buf_2, 2 * w, n, k1, k2 - k1, c1, c2, temp_buf);
        let r = &mut *temp_buf;
        for (i, (a, b)) in std::iter::zip(buf_1.chunks_mut(n), buf_2.chunks_mut(n)).enumerate() {
            let (b1, b2) = b.split_at(w * i);
            let (r1, r2) = r.split_at_mut(n - w * i);
            r1.copy_from_slice(b2);
            for (x, &y) in std::iter::zip(r2.iter_mut(), b1.iter()) {
                *x = -y;
            }
            for ((&x, &y), z) in std::iter::zip(std::iter::zip(a.iter(), r.iter()), b.iter_mut()) {
                *z = c2 * (x - y);
            }
            for (x, &y) in std::iter::zip(a.iter_mut(), r.iter()) {
                *x = c2 * (*x + y);
            }
        }
    }
}


fn ssa_mul<T>(a: &[T], b: &[T]) -> Vec<T>
where
    T: Coefficient,
{
    let r_len = a.len() + b.len() - 1;
    debug_assert!(a.len() > 0 && b.len() > 0 && r_len > 1);
    let log_r_len = ((r_len - 1).ilog2() as usize) + 1;
    let k_ = log_r_len / 2 + 1;
    let d_ = log_r_len - k_;
    
    let n_ = d_ + 1;
    let k = 1 << k_;
    let d = 1 << d_;
    let n = 1 << n_;
    let mut a2 = vec![T::zero(); k * n];
    let mut b2 = vec![T::zero(); k * n];
    
    let w_ = n_ + 1 - k_;
    let w = 1 << w_;
    for (x, y) in std::iter::zip(a.chunks(d), a2.chunks_mut(n)) {
        y[..x.len()].copy_from_slice(x);
    }
    for (x, y) in std::iter::zip(b.chunks(d), b2.chunks_mut(n)) {
        y[..x.len()].copy_from_slice(x);
    }
    
    let ka = a.len().div_ceil(d);
    let kb = b.len().div_ceil(d);
    let kr = ka + kb - 1;
    truncated_ntt(&mut a2, w, n, k, kr, ka);
    truncated_ntt(&mut b2, w, n, k, kr, kb);
    
    if n < SSA_THRESHOLD {
        let mut a3 = vec![T::zero(); n];
        for (r3, b3) in std::iter::zip(a2.chunks_mut(n), b2.chunks(n)).take(kr) {
            a3.copy_from_slice(r3);
            r3.fill(T::zero());
            for (i, &x) in a3.iter().enumerate() {
                let (b4, b5) = b3.split_at(n - i);
                let (r4, r5) = r3.split_at_mut(i);
                for (&y, z) in std::iter::zip(b4.iter(), r5.iter_mut()) {
                    *z += x * y;
                }
                for (&y, z) in std::iter::zip(b5.iter(), r4.iter_mut()) {
                    *z -= x * y;
                }
            }
        }
    }
    else {
        let k2 = n_.div_ceil(2);
        let d2 = n_ - k2;
        for (a3, b3) in std::iter::zip(a2.chunks_mut(n), b2.chunks(n)).take(kr) {
            let r3 = ssa_mul_recursion(a3, b3, k2, d2);
            a3.copy_from_slice(&r3);
        }
    }
    
    let mut r2 = a2;
    truncated_intt(&mut r2, w, n, k, kr);
    
    let mut r = vec![T::zero(); r_len];
    for (i, x) in r2.chunks(n).take(kr).enumerate() {
        let y = &mut r[d * i..std::cmp::min(d * i + n, r_len)];
        for (&x_, y_) in std::iter::zip(x.iter(), y.iter_mut()) {
            *y_ += x_;
        }
    }
    
    r
}

fn ssa_mul_recursion<T>(a: &[T], b: &[T], k_: usize, d_: usize) -> Vec<T>
where
    T: Coefficient,
{
    let n_ = d_ + 1;
    debug_assert!(k_ <= n_ && k_ >= n_ - 1 && k_ >= 1);
    let k = 1 << k_;
    let d = 1 << d_;
    let n = 1 << n_;
    let mut a2 = vec![T::zero(); k * n];
    let mut b2 = vec![T::zero(); k * n];
    
    let u_ = n_ - k_;
    let w_ = u_ + 1;
    let u = 1 << u_;
    let w = 1 << w_;
    for (i, (x, y)) in std::iter::zip(a.chunks(d), a2.chunks_mut(n)).enumerate() {
        if u * i + x.len() <= y.len() {
            y[u * i..u * i + x.len()].copy_from_slice(x);
        }
        else {
            let (x1, x2) = x.split_at(y.len() - u * i);
            y[u * i..].copy_from_slice(x1);
            for (&x_, y_) in std::iter::zip(x2.iter(), y.iter_mut()) {
                *y_ = -x_;
            }
        }
    }
    for (i, (x, y)) in std::iter::zip(b.chunks(d), b2.chunks_mut(n)).enumerate() {
        if u * i + x.len() <= y.len() {
            y[u * i..u * i + x.len()].copy_from_slice(x);
        }
        else {
            let (x1, x2) = x.split_at(y.len() - u * i);
            y[u * i..].copy_from_slice(x1);
            for (&x_, y_) in std::iter::zip(x2.iter(), y.iter_mut()) {
                *y_ = -x_;
            }
        }
    }
    
    ntt(&mut a2, w, n, k);
    ntt(&mut b2, w, n, k);
    
    if n < SSA_THRESHOLD {
        let mut a3 = vec![T::zero(); n];
        for (r3, b3) in std::iter::zip(a2.chunks_mut(n), b2.chunks(n)) {
            a3.copy_from_slice(r3);
            r3.fill(T::zero());
            for (i, &x) in a3.iter().enumerate() {
                let (b4, b5) = b3.split_at(n - i);
                let (r4, r5) = r3.split_at_mut(i);
                for (&y, z) in std::iter::zip(b4.iter(), r5.iter_mut()) {
                    *z += x * y;
                }
                for (&y, z) in std::iter::zip(b5.iter(), r4.iter_mut()) {
                    *z -= x * y;
                }
            }
        }
    }
    else {
        let k2 = n_.div_ceil(2);
        let d2 = n_ - k2;
        for (a3, b3) in std::iter::zip(a2.chunks_mut(n), b2.chunks(n)) {
            let r3 = ssa_mul_recursion(a3, b3, k2, d2);
            a3.copy_from_slice(&r3);
        }
    }
    
    let mut r2 = a2;
    intt(&mut r2, w, n, k);
    
    let mut r = vec![T::zero(); k * d];
    for (i, x) in r2.chunks(n).take(k - 1).enumerate() {
        let y = &mut r[d * i..d * i + n];
        let t = u * i;
        let (x1, x2) = x.split_at(t);
        let (y2, y1) = y.split_at_mut(n - t);
        for (&x_, y_) in std::iter::zip(x1.iter(), y1.iter_mut()) {
            *y_ -= x_;
        }
        for (&x_, y_) in std::iter::zip(x2.iter(), y2.iter_mut()) {
            *y_ += x_;
        }
    }
    let t = u * (k - 1);
    let (_, x) = r2.split_at((k - 1) * n);
    let (x1, x2) = x.split_at(t);
    let (x3, x4) = x1.split_at(t - d);
    let (y1, y2) = r.split_at_mut((k - 1) * d);
    let (y3, _) = y1.split_at_mut(d);
    let (y4, y5) = y2.split_at_mut(n - t);
    for (&x_, y_) in std::iter::zip(x3.iter(), y5.iter_mut()) {
        *y_ -= x_;
    }
    for (&x_, y_) in std::iter::zip(x4.iter(), y3.iter_mut()) {
        *y_ += x_;
    }
    for (&x_, y_) in std::iter::zip(x2.iter(), y4.iter_mut()) {
        *y_ += x_;
    }
    
    r
}


fn ssa_square<T>(a: &[T]) -> Vec<T>
where
    T: Coefficient,
{
    let r_len = 2 * a.len() - 1;
    debug_assert!(a.len() > 0 && r_len > 1);
    let log_r_len = ((r_len - 1).ilog2() as usize) + 1;
    let k_ = log_r_len / 2 + 1;
    let d_ = log_r_len - k_;
    
    let n_ = d_ + 1;
    let k = 1 << k_;
    let d = 1 << d_;
    let n = 1 << n_;
    let mut a2 = vec![T::zero(); k * n];
    
    let w_ = n_ + 1 - k_;
    let w = 1 << w_;
    for (x, y) in std::iter::zip(a.chunks(d), a2.chunks_mut(n)) {
        y[..x.len()].copy_from_slice(x);
    }
    
    let ka = a.len().div_ceil(d);
    let kr = 2 * ka - 1;
    truncated_ntt(&mut a2, w, n, k, kr, ka);
    
    if n < SSA_THRESHOLD {
        let mut a3 = vec![T::zero(); n];
        for r3 in a2.chunks_mut(n).take(kr) {
            a3.copy_from_slice(r3);
            r3.fill(T::zero());
            for (i, &x) in a3.iter().enumerate() {
                let (a4, a5) = a3.split_at(n - i);
                let (r4, r5) = r3.split_at_mut(i);
                for (&y, z) in std::iter::zip(a4.iter(), r5.iter_mut()) {
                    *z += x * y;
                }
                for (&y, z) in std::iter::zip(a5.iter(), r4.iter_mut()) {
                    *z -= x * y;
                }
            }
        }
    }
    else {
        let k2 = n_.div_ceil(2);
        let d2 = n_ - k2;
        for a3 in a2.chunks_mut(n).take(kr) {
            let r3 = ssa_square_recursion(a3, k2, d2);
            a3.copy_from_slice(&r3);
        }
    }
    
    let mut r2 = a2;
    truncated_intt(&mut r2, w, n, k, kr);
    
    let mut r = vec![T::zero(); r_len];
    for (i, x) in r2.chunks(n).take(kr).enumerate() {
        let y = &mut r[d * i..std::cmp::min(d * i + n, r_len)];
        for (&x_, y_) in std::iter::zip(x.iter(), y.iter_mut()) {
            *y_ += x_;
        }
    }
    
    r
}

fn ssa_square_recursion<T>(a: &[T], k_: usize, d_: usize) -> Vec<T>
where
    T: Coefficient,
{
    let n_ = d_ + 1;
    debug_assert!(k_ <= n_ && k_ >= n_ - 1 && k_ >= 1);
    let k = 1 << k_;
    let d = 1 << d_;
    let n = 1 << n_;
    let mut a2 = vec![T::zero(); k * n];
    
    let u_ = n_ - k_;
    let w_ = u_ + 1;
    let u = 1 << u_;
    let w = 1 << w_;
    for (i, (x, y)) in std::iter::zip(a.chunks(d), a2.chunks_mut(n)).enumerate() {
        if u * i + x.len() <= y.len() {
            y[u * i..u * i + x.len()].copy_from_slice(x);
        }
        else {
            let (x1, x2) = x.split_at(y.len() - u * i);
            y[u * i..].copy_from_slice(x1);
            for (&x_, y_) in std::iter::zip(x2.iter(), y.iter_mut()) {
                *y_ = -x_;
            }
        }
    }
    
    ntt(&mut a2, w, n, k);
    
    if n < SSA_THRESHOLD {
        let mut a3 = vec![T::zero(); n];
        for r3 in a2.chunks_mut(n) {
            a3.copy_from_slice(r3);
            r3.fill(T::zero());
            for (i, &x) in a3.iter().enumerate() {
                let (a4, a5) = a3.split_at(n - i);
                let (r4, r5) = r3.split_at_mut(i);
                for (&y, z) in std::iter::zip(a4.iter(), r5.iter_mut()) {
                    *z += x * y;
                }
                for (&y, z) in std::iter::zip(a5.iter(), r4.iter_mut()) {
                    *z -= x * y;
                }
            }
        }
    }
    else {
        let k2 = n_.div_ceil(2);
        let d2 = n_ - k2;
        for a3 in a2.chunks_mut(n) {
            let r3 = ssa_square_recursion(a3, k2, d2);
            a3.copy_from_slice(&r3);
        }
    }
    
    let mut r2 = a2;
    intt(&mut r2, w, n, k);
    
    let mut r = vec![T::zero(); k * d];
    for (i, x) in r2.chunks(n).take(k - 1).enumerate() {
        let y = &mut r[d * i..d * i + n];
        let t = u * i;
        let (x1, x2) = x.split_at(t);
        let (y2, y1) = y.split_at_mut(n - t);
        for (&x_, y_) in std::iter::zip(x1.iter(), y1.iter_mut()) {
            *y_ -= x_;
        }
        for (&x_, y_) in std::iter::zip(x2.iter(), y2.iter_mut()) {
            *y_ += x_;
        }
    }
    let t = u * (k - 1);
    let (_, x) = r2.split_at((k - 1) * n);
    let (x1, x2) = x.split_at(t);
    let (x3, x4) = x1.split_at(t - d);
    let (y1, y2) = r.split_at_mut((k - 1) * d);
    let (y3, _) = y1.split_at_mut(d);
    let (y4, y5) = y2.split_at_mut(n - t);
    for (&x_, y_) in std::iter::zip(x3.iter(), y5.iter_mut()) {
        *y_ -= x_;
    }
    for (&x_, y_) in std::iter::zip(x4.iter(), y3.iter_mut()) {
        *y_ += x_;
    }
    for (&x_, y_) in std::iter::zip(x2.iter(), y4.iter_mut()) {
        *y_ += x_;
    }
    
    r
}


#[derive(Clone, PartialEq, Eq, Debug)]
struct SparsePolynomial<T>(BTreeMap<usize, T>);

impl<T> SparsePolynomial<T>
where
    T: Coefficient,
{
    pub fn degree(&self) -> Option<usize> {
        self.0.last_key_value().map(|(&d, _)| d)
    }
    
    pub fn monic(&self) -> SparsePolynomial<T> {
        let mut a = self.clone();
        
        let Some((_, &l)) = a.0.last_key_value() else {
            return a;
        };
        
        let c = l.inv();
        for x in a.0.values_mut() {
            *x *= c;
        }
        
        a
    }
}

impl<T> From<BTreeMap<usize, T>> for SparsePolynomial<T>
where
    T: Coefficient,
{
    fn from(x: BTreeMap<usize, T>) -> SparsePolynomial<T> {
        let mut x = x;
        x.retain(|_, &mut c| !c.is_zero());
        SparsePolynomial(x)
    }
}

impl<T> From<SparsePolynomial<T>> for Polynomial<T>
where
    T: Coefficient,
{
    fn from(x: SparsePolynomial<T>) -> Polynomial<T> {
        let n = x.degree().map_or(0, |x| x + 1);
        let mut v = vec![T::zero(); n];
        for (&e, &c) in x.0.iter() {
            v[e] = c;
        }
        Polynomial(v)
    }
}

impl<T> Add<T> for &SparsePolynomial<T>
where
    T: Coefficient,
{
    type Output = SparsePolynomial<T>;
    
    fn add(self, other: T) -> SparsePolynomial<T> {
        let (mut a, b) = (self.0.clone(), other);
        let x = a.entry(0).or_insert(T::zero());
        *x += b;
        SparsePolynomial::from(a)
    }
}

impl<T> Sub<T> for &SparsePolynomial<T>
where
    T: Coefficient,
{
    type Output = SparsePolynomial<T>;
    
    fn sub(self, other: T) -> SparsePolynomial<T> {
        let (mut a, b) = (self.0.clone(), other);
        let x = a.entry(0).or_insert(T::zero());
        *x -= b;
        SparsePolynomial::from(a)
    }
}

impl<T> Rem<&SparsePolynomial<T>> for Polynomial<T>
where
    T: Coefficient,
{
    type Output = Polynomial<T>;
    
    fn rem(self, other: &SparsePolynomial<T>) -> Polynomial<T> {
        let (a, b) = (self, other);
        let Some(d) = b.degree() else {
            panic!();
        };
        
        let Some(n) = a.degree().and_then(|x| x.checked_sub(d)) else {
            return a;
        };
        
        let b = b.monic();
        let mut b2 = b.0;
        b2.pop_last();
        
        let mut a2 = a.0;
        for i in (0..=n).rev() {
            let c = a2[d + i];
            for (&k, &v) in b2.iter() {
                a2[k + i] -= v * c;
            }
        }
        a2.truncate(d);
        
        Polynomial::from(a2)
    }
}

impl<T> RemAssign<&SparsePolynomial<T>> for Polynomial<T>
where
    T: Coefficient,
{
    fn rem_assign(&mut self, other: &SparsePolynomial<T>) {
        let v = std::mem::take(&mut self.0);
        let (a, b) = (Polynomial(v), other);
        *self = a.rem(b);
    }
}


fn main() {
    use std::time::Instant;
    
    const P: u64 = 59_u64.wrapping_neg();
    
    let t0 = Instant::now();
    
    #[allow(non_snake_case)]
    let E = [(1 << 24) + 17, (1 << 24) + 3, 3, 2, 1, 0];
    #[allow(non_snake_case)]
    let C = [1, 17761542461647558231, 13293668011354679701, 9204760597720472707, 8540722934676348527, 3568330912555059249];
    let pash = b"\x29\x77\x9e\x08\x3d\x3a\x77\xa1\xa7\x8a\x1c\xe1\x45\xc0\x89\x0c\x59\x96\xbd\x91\xd1\x1c\xdf\x9f\x04\x9d\xe5\x89\x5f\xae\x1d\xe6\x12\xf2\x37\xcc\x64\x88\xea\xb2\x79\x65\xd4\xcb\x2e\x31\x39\xaf\xfb\x83\xc7\x90\x35\x26\x5e\x08\xb5\x0e\xed\x5e\xed\x8a\xd3\x74\x1e\xca\x2a\x42\x81\x40\x57\xe9\x2e\xc5\xd8\x57\xc6\xc1\x7e\x4a\x95\xe7\xfe\xc3\x9e\x6d\xea\x49\xa9\xe3\x69\x35\x4b\xa9\xf7\x2c\x13\x8d\x5e\x0a\x5a\x81\x37\x41\x08\x0e\x87\x98\x23\x82\x93\x2f\xfe\x51\xe0\x25\xad\x4d\xb7\x11\x47\xa6\xcf\xde\x90\x82\x1d\x14\x3c\xb5\xe8\x8f\x57\x42\x83\x4b\x2a\x52\x7c\xb6\x31\x3e\x92\x81\x40\x4d\x5f\xfb\xd5\xbe\xcd\x54\xb2\xf7\x4a\xeb\x89\xf7\x07\x9c\x23\xcf\x9c\x9e\x2b\x6d\xa2\x89\x63\x13\xae\xef\x92\x82\x56\x7f\xad\xb3\x10\xb6\x1f\x30\x9e\xfb\xfa\xb5\xd6\x2c\xbe\x65\x0c\x81";
    
    let pash: Vec<_> = pash.chunks_exact(8).map(|x| u64::from_be_bytes(x.try_into().unwrap())).collect();
    let f0 = SparsePolynomial::from(std::iter::zip(E.iter(), C.iter()).map(|(&e, &c)| (e, MontgomeryForm::<P>::from(c))).collect::<BTreeMap<_, _>>());
    
    let mut flag = vec![];
    for (idx, &h) in pash.iter().enumerate() {
        let t1 = Instant::now();
        let f1 = &f0 - MontgomeryForm::<P>::from(h);
        
        let d = P as usize;
        let n = d.checked_ilog2().map_or(0, |x| x + 1);
        let l = 25;
        let n2 = n.saturating_sub(l);
        let d2 = d / (1 << n2);
        let mut f2 = (Polynomial::one() << d2) % &f1;
        for i in (0..n2).rev() {
            f2 = f2.square();
            if (d >> i) & 0x1 != 0 {
                f2 <<= 1;
            }
            f2 %= &f1;
        }
        f2 -= Polynomial::x();
        
        let f1 = Polynomial::from(f1);
        let f3 = f2.gcd(&f1).monic();
        
        let mut factors = vec![(f3, 0)];
        let mut roots = vec![];
        let h = (P - 1) / 2;
        let n = h.checked_ilog2().map_or(0, |x| x + 1);
        while let Some((f4, d)) = factors.pop() {
            if !f4.degree().is_some_and(|x| x > 1) {
                if f4.degree().is_some_and(|x| x == 1) {
                    roots.push(-f4.coefficient(0));
                }
                continue;
            }
            
            let mut f5 = Polynomial::one();
            let f6 = Polynomial::from([d, 1].iter().map(|&x| MontgomeryForm::<P>::from(x)).collect::<Vec<_>>());
            for i in (0..n).rev() {
                f5 = f5.square();
                if (h >> i) & 0x1 != 0 {
                    f5 *= &f6;
                }
                f5 %= &f4;
            }
            f5 -= Polynomial::one();
            
            let f6 = f5.gcd(&f4).monic();
            factors.push((f4 / &f6, d + 1));
            factors.push((f6, d + 1));
        }
        
        let t = t1.elapsed();
        let roots: Vec<_> = roots.iter().map(|&x| u64::from(x)).collect();
        let pt: Vec<_> = roots.iter().map(|&x| x.to_be_bytes()).collect();
        println!("{:2}/{:2}: t = {:.1?}, pt = [{}]", idx + 1, pash.len(), t, pt.iter().map(|&x| format!("b\"{}\"", x.escape_ascii())).collect::<Vec<_>>().join(", "));
        if let Some(b) = pt.iter().max_by_key(|&b| b.iter().filter(|&&x| x.is_ascii()).count()) {
            flag.extend(b);
        }
    }
    
    println!("b\"{}\"", flag.escape_ascii());
    
    let t = t0.elapsed();
    println!("{:?}", t);
}
