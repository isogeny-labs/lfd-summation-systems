//! Open-addressing hash map u64 -> u64 (keys < 2^63; u64::MAX = empty).
pub struct HMap {
    keys: Vec<u64>,
    vals: Vec<u64>,
    shift: u32,
    pub len: usize,
}
const EMPTY: u64 = u64::MAX;
impl HMap {
    pub fn with_capacity(n: usize) -> Self {
        let cap = (n.max(4) * 2).next_power_of_two();
        HMap { keys: vec![EMPTY; cap], vals: vec![0; cap], shift: 64 - cap.trailing_zeros(), len: 0 }
    }
    #[inline(always)]
    fn slot(&self, k: u64) -> usize { (k.wrapping_mul(0x9E3779B97F4A7C15) >> self.shift) as usize }
    /// insert if absent; returns true if newly inserted
    pub fn insert(&mut self, k: u64, v: u64) -> bool {
        if (self.len + 1) * 2 > self.keys.len() { self.grow(); }
        let mask = self.keys.len() - 1;
        let mut i = self.slot(k);
        loop {
            let kk = self.keys[i];
            if kk == EMPTY { self.keys[i] = k; self.vals[i] = v; self.len += 1; return true; }
            if kk == k { return false; }
            i = (i + 1) & mask;
        }
    }
    #[inline(always)]
    pub fn get(&self, k: u64) -> Option<u64> {
        let mask = self.keys.len() - 1;
        let mut i = self.slot(k);
        loop {
            let kk = unsafe { *self.keys.get_unchecked(i) };
            if kk == k { return Some(unsafe { *self.vals.get_unchecked(i) }); }
            if kk == EMPTY { return None; }
            i = (i + 1) & mask;
        }
    }
    #[inline(always)]
    pub fn contains(&self, k: u64) -> bool { self.get(k).is_some() }
    fn grow(&mut self) {
        let old_k = std::mem::take(&mut self.keys);
        let old_v = std::mem::take(&mut self.vals);
        let cap = old_k.len() * 2;
        self.keys = vec![EMPTY; cap];
        self.vals = vec![0; cap];
        self.shift = 64 - cap.trailing_zeros();
        self.len = 0;
        for (k, v) in old_k.into_iter().zip(old_v) { if k != EMPTY { self.insert(k, v); } }
    }
    pub fn bytes(&self) -> usize { self.keys.len() * 16 }
}
