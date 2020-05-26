extern crate minifb;

use minifb::{Key, Window, WindowOptions, MouseMode, MouseButton};
use std::time::Instant;

const WIDTH: usize = 640;
const HEIGHT: usize = 480;
const MAX_RAY_DEPTH: u32 = 5;


#[derive(Debug,Copy,Clone)]
struct Vec3 {
    x: f32,
    y: f32,
    z: f32,
}

impl Vec3 {
    fn new(x: f32, y: f32, z:f32 ) -> Vec3 {
        Vec3 { x: x, y: y, z: z }
    }

    fn new_const(x: f32) -> Vec3 {
        Vec3 { x: x, y: x, z: x }
    }

    fn dot(&self, other: Vec3) -> f32
    {
        self.x * other.x +
        self.y * other.y +
        self.z * other.z
    }

    fn length2(&self) -> f32
    {
        self.x * self.x +
        self.y * self.y +
        self.z * self.z
    }

    fn normalize(&mut self)
    {
        let norm = self.length2().sqrt();
        self.x /= norm;
        self.y /= norm;
        self.z /= norm;
        return
    }

    fn to_color(&self) -> u32 {
        //assert!(self.x <= 1.0);
        //assert!(self.y <= 1.0);
        //assert!(self.z <= 1.0);
        (((self.x * 255.99) as u32) << 16) |
        (((self.y * 255.99) as u32) <<  8) |
         ((self.z * 255.99) as u32)
    }
}

impl std::ops::Add for Vec3
{
    type Output = Vec3;
    fn add(self, rhs: Vec3) -> Self::Output {
        Vec3 { x: self.x + rhs.x,
               y: self.y + rhs.y,
               z: self.z + rhs.z }
    }
}

impl std::ops::Mul for Vec3
{
    type Output = Vec3;
    fn mul(self, rhs: Vec3) -> Self::Output {
        Vec3 { x: self.x * rhs.x,
               y: self.y * rhs.y,
               z: self.z * rhs.z }
    }
}

impl std::ops::Mul<f32> for Vec3
{
    type Output = Vec3;
    fn mul(self, rhs: f32) -> Self::Output {
        Vec3 { x: self.x * rhs,
               y: self.y * rhs,
               z: self.z * rhs }
    }
}

impl std::ops::Sub for Vec3
{
    type Output = Vec3;
    fn sub(self, rhs: Vec3) -> Self::Output {
        Vec3 { x: self.x - rhs.x,
               y: self.y - rhs.y,
               z: self.z - rhs.z }
    }
}

impl std::ops::Neg for Vec3
{
    type Output = Vec3;
    fn neg(self) -> Self::Output {
        Vec3 { x: -self.x,
               y: -self.y,
               z: -self.z }
    }
}

impl std::ops::AddAssign for Vec3 {
    fn add_assign(&mut self, rhs: Vec3) {
        *self = Self { x: self.x + rhs.x,
                       y: self.y + rhs.y,
                       z: self.z + rhs.z };
    }
}

struct Sphere {
    center: Vec3,
    radius2: f32,
    surface_color: Vec3,
    emission_color: Vec3,
    transparency: f32,
    reflection: f32,
}

impl Sphere {
    fn new(center: Vec3, radius: f32, sc: Vec3,
           reflect: f32, trans: f32) -> Sphere {
        let ec = Vec3::new_const(0.0);
        Sphere::new_withlight(center, radius, sc, reflect, trans, ec)
    }

    fn new_withlight(center: Vec3, radius: f32, sc: Vec3,
           reflect: f32, trans: f32, ec: Vec3) -> Sphere {
        Sphere { center: center, radius2: radius*radius,
                 surface_color: sc, reflection: reflect,
                 transparency: trans, emission_color: ec }
    }

    fn intersect(&self, rayorig: Vec3, raydir: Vec3) -> (bool, f32, f32) {
        let l = self.center - rayorig;
        let tca = l.dot(raydir);
        if tca < 0.0 {
            return (false, 0.0, 0.0);
        }
        let d2 = l.dot(l) - tca*tca;
        if d2 > self.radius2 {
            return (false, 0.0, 0.0);
        }
        let thc = (self.radius2 - d2).sqrt();
        let t0 = tca - thc;
        let t1 = tca + thc;
        (true, t0, t1)
    }
}

fn mix(a: f32, b: f32, m: f32) -> f32
{
    b * m + a * (1.0 - m)
}

fn trace(rayorig: Vec3, raydir: Vec3, spheres: &Vec<Sphere>, depth: u32, debug: bool) -> Vec3
{
    let mut tnear = f32::INFINITY;
    let mut maybe_sphere : Option<&Sphere> = None;
    for s in spheres {
        let (found, mut t0, t1) = s.intersect(rayorig, raydir);
        if found {
            if t0 < 0.0 { t0 = t1; }
            if t0 < tnear {
                tnear = t0;
                maybe_sphere = Some(&s);
            }
        }
    }
    if maybe_sphere.is_none() {
        return Vec3::new_const(2.0);
    }
    let sphere = maybe_sphere.unwrap();

    let mut surface_color = Vec3::new_const(0.0);
    let phit = rayorig + raydir * tnear;
    let mut nhit = phit - sphere.center;
    nhit.normalize();

    if debug {
        println!("nhit: {:?}", nhit);
    }

    let mut inside = false;
    let bias = 1e-4_f32;
    if raydir.dot(nhit) > 0.0 {
        nhit = -nhit;
        inside = true;
    }

    if (sphere.transparency > 0.0 || sphere.reflection > 0.0) && depth < MAX_RAY_DEPTH {
        let facingratio = -raydir.dot(nhit);
        let fresneleffect = mix((1.0-facingratio).powf(3.0), 1.0, 0.1);
        let mut refldir = raydir - nhit * 2.0 * raydir.dot(nhit);
        refldir.normalize();
        if debug {
            println!("refldir: {:?}", refldir);
            println!("facingratio: {}", facingratio);
            println!("fresneleffect: {}", fresneleffect);
        }
        let reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1, debug);
        let mut refraction = Vec3::new_const(0.0);
        if sphere.transparency != 0.0 {
            let ior = 1.1;
            let eta = if inside { ior } else { 1.0 / ior };
            let cosi = -nhit.dot(raydir);
            let k = 1.0 - eta * eta * (1.0 - cosi * cosi);
            let mut refrdir = raydir * eta + nhit * (eta * cosi - k.sqrt());
            refrdir.normalize();
            refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1, debug);
            if debug {
                println!("refrdir: {:?}", refrdir);
                println!("refraction: {:?}", refraction);
            }
        }
        surface_color = (
            reflection * fresneleffect +
            refraction * (1.0 - fresneleffect) * sphere.transparency) * sphere.surface_color;
        if debug {
            println!("surface_color: {:?}", surface_color);
            println!("reflection: {:?} fresneleffect: {:?} \
                     refraction: {:?} sphere.transparency: {:?} \
                     sphere.surface_color: {:?}",
                     reflection, fresneleffect, refraction,
                     sphere.transparency, sphere.surface_color);
        }
    }
    else {
        for (i, s1) in spheres.iter().enumerate() {
            if s1.emission_color.x > 0.0 {
                // this is a light
                let mut transmission = Vec3::new_const(1.0);
                let mut light_direction = s1.center - phit;
                light_direction.normalize();
                for (j,s2) in spheres.iter().enumerate() {
                    if i != j {
                        let (hit, _, _) = s2.intersect(phit + nhit * bias, light_direction);
                        if hit {
                            transmission = Vec3::new_const(0.0);
                            break;
                        }
                    }
                }
                let max = 0.0_f32.max(nhit.dot(light_direction));
                surface_color += sphere.surface_color * transmission * max * s1.emission_color;
            }
        }
    }
    surface_color + sphere.emission_color
}

fn render(spheres: &Vec<Sphere>, buffer: &mut Vec<u32>) {
    let invwidth = 1.0 / (WIDTH as f32);
    let invheight = 1.0 / (HEIGHT as f32);
    let fov = 30.0;
    let aspectratio = (WIDTH as f32)/ (HEIGHT as f32);
    let angle = (std::f32::consts::PI * 0.5 * fov / 180.0).tan();
    for y in 0..HEIGHT {
        for x in 0..WIDTH {
            let xx = (2.0 * ((x as f32 + 0.5) * invwidth) - 1.0) * angle * aspectratio;
            let yy = (1.0 - 2.0 * ((y as f32 + 0.5) * invheight)) * angle;
            let mut raydir = Vec3::new(xx, yy, -1.0);
            raydir.normalize();
            let traced = trace(Vec3::new_const(0.0), raydir, spheres, 0, x == 320 && y == 366 && false);
            if x == 320 && y == 366 && false {
                println!("xx: {} yy: {} raydir: {:?}", xx, yy, raydir);
                println!("Debug: 320,366 -> {:?}", traced);
            }
            buffer[y*WIDTH+x] = traced.to_color();
        }
    }
}

fn main() {
    let mut buffer: Vec<u32> = vec![0; WIDTH * HEIGHT];

    let mut window = Window::new(
        "RayBan - ESC to exit",
        WIDTH,
        HEIGHT,
        WindowOptions::default(),
    )
    .unwrap_or_else(|e| {
        panic!("{}", e);
    });

    let spheres = vec![
        Sphere::new(Vec3::new(0.0,-10004.0,-20.0),  10000.0, Vec3::new(0.20,0.20,0.20), 0.0, 0.0),
        Sphere::new(Vec3::new(0.0,0.0,-20.0),           4.0, Vec3::new(1.00,0.32,0.36), 1.0, 0.5),
        Sphere::new(Vec3::new(5.0,-1.0,-15.0),          2.0, Vec3::new(0.90,0.76,0.46), 1.0, 0.0),
        Sphere::new(Vec3::new(5.0,0.0,-25.0),           3.0, Vec3::new(0.65,0.77,0.97), 1.0, 0.0),
        Sphere::new(Vec3::new(-5.5,0.0,-15.0),          3.0, Vec3::new(0.90,0.90,0.90), 1.0, 0.0),
        // light
        Sphere::new_withlight(Vec3::new(0.0,20.0,-30.0), 3.0, Vec3::new(0.00,0.00,0.00 ), 0.0, 0.0,
            Vec3::new_const(3.0)),
    ];

    // Limit to max ~60 fps update rate
    window.limit_update_rate(Some(std::time::Duration::from_micros(16600)));

    while window.is_open() && !window.is_key_down(Key::Escape) {
        let start = Instant::now();
        render(&spheres, &mut buffer);
        let ms = start.elapsed().as_millis();
        println!("Rendered in {:?} ms", ms);

        if window.get_mouse_down(MouseButton::Left) {
            let (x,y) = window.get_mouse_pos(MouseMode::Clamp).unwrap();
            let xx = x as usize;
            let yy = y as usize;
            //println!("{}, {} ({},{})", x, y, xx, yy);
            println!("({},{}) = {:#x}", xx, yy, buffer[yy*WIDTH+xx]);
            break;
        }

        // We unwrap here as we want this code to exit if it fails. Real applications may want to handle this in a different way
        window
            .update_with_buffer(&buffer, WIDTH, HEIGHT)
            .unwrap();
    }
}
