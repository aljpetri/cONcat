//Interval represents an interval of a read that is covered by a fragment
#[derive(Debug, PartialEq,Eq,Clone,PartialOrd,Ord)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
}
/*impl Ord for Interval {
    fn cmp(&self, other: &Self) -> Ordering {
        self.start. cmp(&other. start)
    }
}

impl PartialEq<Self> for Interval {
    fn eq(&self, other: &Self) -> bool {
        if self.start == other.start && self.end == other.end{
            return true
        }
        else {
            return false
        }
    }
}*/


