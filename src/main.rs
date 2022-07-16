use std::collections::HashMap;

fn main() {
    let iterations: u64 = 8;

    let vars: Vars = Vars {
        horiz_step: iterations,
        vert_step: iterations + 2,
        s_min: 0.0,
        s_max: 40.0,
        k: 21.0,
        t: 0.33333,
        r: 0.1,
        v: 0.4,
        q: 0.0
    };
    
    let abc: ABC = ABC {
        arr_len: vars.vert_step,
        dt: vars.calc_dt(),
        r: vars.r,
        v: vars.v
    };
    
    let mut ans: FDEngine = FDEngine {
        horiz_step: vars.horiz_step,
        ds_list: vars.calc_ds_list(),
        abc: abc.get_array(),
        top_layer: vars.calc_top_boundary(),
        counter: 0
    };

    let bound: Vec<f64> = vars.calc_end_boundary();
    
    println!("{:#?}", ans.evaluate(bound));
}

pub struct Vars {
    horiz_step: u64,
    vert_step: u64,
    s_min: f64,
    s_max: f64,
    k: f64,
    t: f64,
    r: f64,
    v: f64,
    q: f64
}

impl Vars {
    pub fn calc_ds(&self) -> f64 {
        // Calculates equidistant finite distance step of PDE
        return (self.s_max - self.s_min) / self.vert_step as f64
    }
    
    pub fn calc_ds_list(&self) -> Vec<f64> {
        // Returns left boundary of lattice in Vector
        let ds: f64 = self.calc_ds();
        let mut boundary: Vec<f64> = Vec::new();
        for n in 0..=self.vert_step {
            boundary.push(n as f64 * ds);
        }
        return boundary
    }

    pub fn calc_dt(&self) -> f64 {
        // Calculates timespan between equidistiant steps
        return self.t / self.horiz_step as f64 
    }

    pub fn calc_end_boundary(&self) -> Vec<f64> {
        // Returns right boundary of lattice in Vector
        let zero: f64 = 0.0;
        let ds: f64 = self.calc_ds();
        let mut boundary: Vec<f64> = Vec::new();
        for bound_res in 0..=self.vert_step {
            boundary.push(zero.max(self.k - bound_res as f64 * ds));
        }
        return boundary 
    }

    pub fn calc_top_boundary(&self) -> Vec<f64> {
        // Returns top boundary of lattice in Vector
        let dt: f64 = self.calc_dt();
        let mut boundary: Vec<f64> = Vec::new();
        for t in 0..self.horiz_step {
            let formula: f64 = 
            self.k * f64::exp(
            - (self.horiz_step as f64 - t as f64) * dt * self.r)
            - (0.0 * f64::exp(
            - (self.horiz_step as f64 - t as f64) * dt * self.q));
            boundary.push(formula);
        }
        return boundary
    }
}

pub struct ABC {
    arr_len: u64,
    dt: f64,
    r: f64,
    v: f64
}

impl ABC {
    pub fn calc_a(&self) -> Vec<f64> {
        // Returns "a vector" of computation variables
        let mut res: Vec<f64> = Vec::new();
        for a in 0..self.arr_len + 1 {
            let formula: f64 = (
            - 0.5 * self.r * self.dt * a as f64
            + 0.5 * f64::powf(self.v, 2.0) * self.dt
            * f64::powf(a as f64, 2.0)) / (1.0 + self.r * self.dt);
            res.push(formula);
        }
        return res
    }
    pub fn calc_b(&self) -> Vec<f64> {
        // Returns "b vector" of computation variables
        let mut res: Vec<f64> = Vec::new();
        for b in 0..self.arr_len + 1 {
            let formula: f64 = 
            1.0 / (1.0 + self.r * self.dt)
            * (1.0 - f64::powf(self.v, 2.0) * f64::powf(b as f64, 2.0)
            * self.dt);
            res.push(formula);
        }
        return res
    }
    pub fn calc_c(&self) -> Vec<f64> {
        // Returns "c vector" of computation variables
        let mut res: Vec<f64> = Vec::new();
        for c in 0..self.arr_len + 1 {
            let formula: f64 = (
            0.5 * self.r * self.dt * c as f64
            + 0.5 * f64::powf(self.v, 2.0) * self.dt
            * f64::powf(c as f64, 2.0)) / (1.0 + self.r * self.dt);
            res.push(formula);
        }
        return res
    }
    pub fn get_array(&self) -> HashMap<&str, Vec<f64>> {
        // Returns ABC computation variables in one run
        let mut res: HashMap<&str, Vec<f64>> = HashMap::new();
        res.insert("a", self.calc_a());
        res.insert("b", self.calc_b());
        res.insert("c", self.calc_c());
        return res
    }
}

pub struct FDEngine<'a> {
    horiz_step: u64,
    ds_list: Vec<f64>,
    abc: HashMap<&'a str, Vec<f64>>,
    top_layer: Vec<f64>,
    counter: u64
}

impl<'a> FDEngine<'a> {
    pub fn compute_mmult(&self, row_arr: Vec<f64>, col_arr: Vec<f64>) -> f64 {
        let mut res: f64 = 0.0;
        for position in 0..row_arr.len() {
            res += row_arr[position] * col_arr[position];
        }
        return res
    }
    pub fn evaluate(&mut self, boundary: Vec<f64>) -> HashMap<String, f64> {
        let mut ans_vec: Vec<f64> = Vec::new();
        if self.counter == self.horiz_step {
            let mut res: HashMap<String, f64> = HashMap::new();
            for value in 0..self.ds_list.len() {
                let index_to_push: String = self.ds_list[value].to_string();
                let value_to_push: f64 = boundary[value];
                res.insert(index_to_push, value_to_push);
            }
            return res
        } else {
            for i in 0..boundary.len() - 1 {
                if i == 0 {
                    ans_vec.push(self.top_layer[self.top_layer.len() - 1]);
                } else {
                    let mmult: f64 = self.compute_mmult(
                    vec![boundary[i - 1], boundary[i], boundary[i + 1]],
                    vec![self.abc["a"][i], self.abc["b"][i], self.abc["c"][i]]
                );
                ans_vec.push(mmult);
                }
            }
        }
        ans_vec.push(0.0);
        self.counter += 1;
        self.top_layer.pop();
        return self.evaluate(ans_vec)
    }
}