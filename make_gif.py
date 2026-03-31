import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import argparse

# 1. Torus Parameters
R = 3.0  # Distance from center of hole to center of tube
r = 1.0  # Radius of the tube
offset = 0.05 # Gap to prevent z-fighting

# Grid for the torus surface
u = np.linspace(0, 2 * np.pi, 40)
v = np.linspace(0, 2 * np.pi, 40)
U, V = np.meshgrid(u, v)

X = (R + (r - offset) * np.cos(V)) * np.cos(U)
Y = (R + (r - offset) * np.cos(V)) * np.sin(U)
Z = (r - offset) * np.sin(V)

# Setup the figure
fig = plt.figure(figsize=(8, 6), facecolor='white')
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor('white')
ax.set_axis_off()

# Draw the torus surface (zorder=1 to allow tracks to easily sit on top)
ax.plot_surface(X, Y, Z, color='#E8E8EE', edgecolors='gray', lw=0.3, alpha=0.3, zorder=1)

# ==========================================
# 2. STATIC TRACKS & TRIVIAL LOOP
# ==========================================

# X_L Track (Red)
u_xl_f = np.linspace(np.pi/2, 1.5 * np.pi, 50) # Front half
X_xl_f = (R + (r + offset) * np.cos(np.pi/2)) * np.cos(u_xl_f)
Y_xl_f = (R + (r + offset) * np.cos(np.pi/2)) * np.sin(u_xl_f)
Z_xl_f = (r + offset) * np.sin(np.pi/2) * np.ones_like(u_xl_f)
# Boosted visibility and locked to front
ax.plot(X_xl_f, Y_xl_f, Z_xl_f, color='#D9534F', lw=5, alpha=1.0, zorder=10)

u_xl_b = np.linspace(-np.pi/2, np.pi/2, 50) # Back half
X_xl_b = (R + (r + offset) * np.cos(np.pi/2)) * np.cos(u_xl_b)
Y_xl_b = (R + (r + offset) * np.cos(np.pi/2)) * np.sin(u_xl_b)
Z_xl_b = (r + offset) * np.sin(np.pi/2) * np.ones_like(u_xl_b)
ax.plot(X_xl_b, Y_xl_b, Z_xl_b, color='#D9534F', lw=4, alpha=1.0, zorder=10)

# Z_L Track (Blue) at Right Side
u_zl = 1.5 * np.pi

# Front Half: Bound precisely from -0.1*pi to 1.1*pi based on new logic
v_zl_f = np.linspace(-0.1 * np.pi, 1.1 * np.pi, 60) 
X_zl_f = (R + (r + offset) * np.cos(v_zl_f)) * np.cos(u_zl)
Y_zl_f = (R + (r + offset) * np.cos(v_zl_f)) * np.sin(u_zl)
Z_zl_f = (r + offset) * np.sin(v_zl_f)
# Boosted visibility and locked to front
ax.plot(X_zl_f, Y_zl_f, Z_zl_f, color='#4A90E2', lw=5, alpha=1.0, zorder=10)

# Back Half: The remaining segment from 1.1*pi to 1.9*pi
v_zl_b = np.linspace(1.1 * np.pi, 1.9 * np.pi, 40)
X_zl_b = (R + (r + offset) * np.cos(v_zl_b)) * np.cos(u_zl)
Y_zl_b = (R + (r + offset) * np.cos(v_zl_b)) * np.sin(u_zl)
Z_zl_b = (r + offset) * np.sin(v_zl_b)
ax.plot(X_zl_b, Y_zl_b, Z_zl_b, color='#4A90E2', lw=4, alpha=0.8, linestyle='--', zorder=0)

# Static Square Trivial Loop
du, dv = 0.5, 0.8 
u0, v0 = np.pi - du/2, np.pi/4 - dv/2
u_sq = np.array([u0, u0+du, u0+du, u0, u0])
v_sq = np.array([v0, v0, v0+dv, v0+dv, v0])

X_triv = (R + (r + offset) * np.cos(v_sq)) * np.cos(u_sq)
Y_triv = (R + (r + offset) * np.cos(v_sq)) * np.sin(u_sq)
Z_triv = (r + offset) * np.sin(v_sq)
ax.plot(X_triv, Y_triv, Z_triv, color='#82B366', lw=5, zorder=10)

ax.set_xlim(-4, 4)
ax.set_ylim(-4, 4)
ax.set_zlim(-3, 3)

# View angle locked
ax.view_init(elev=35., azim=180)

# ==========================================
# 3. ANIMATION LOGIC (MANUAL ARROWS)
# ==========================================
arrow_artists = []
total_frames = 60

def update(frame):
    global arrow_artists
    for artist in arrow_artists:
        artist.remove()
    arrow_artists.clear()
        
    prog = (frame / total_frames) * 2 * np.pi
    
    # ------------------------------------------------
    # Red Arrow (X_L)
    # ------------------------------------------------
    u_r = prog
    v_r = np.pi / 2
    
    p_r = np.array([(R + (r + offset) * np.cos(v_r)) * np.cos(u_r),
                    (R + (r + offset) * np.cos(v_r)) * np.sin(u_r),
                    (r + offset) * np.sin(v_r)])
    
    v_vec_r = np.array([-np.sin(u_r), np.cos(u_r), 0])
    w_vec_r = np.array([-np.cos(u_r), -np.sin(u_r), 0])
    
    tip_r = p_r + v_vec_r * 1.0
    base_r = tip_r - v_vec_r * 0.4
    left_r = base_r + w_vec_r * 0.25
    right_r = base_r - w_vec_r * 0.25
    
    # RED ARROW: Always on top as requested
    red_layer = 100 
    
    l1, = ax.plot([p_r[0], tip_r[0]], [p_r[1], tip_r[1]], [p_r[2], tip_r[2]], color="#942723", lw=4, zorder=red_layer)
    l2, = ax.plot([left_r[0], tip_r[0], right_r[0]], [left_r[1], tip_r[1], right_r[1]], [left_r[2], tip_r[2], right_r[2]], color="#942723", lw=4, zorder=red_layer)
    arrow_artists.extend([l1, l2])

    # ------------------------------------------------
    # Blue Arrow (Z_L)
    # ------------------------------------------------
    u_b = 1.5 * np.pi
    v_b = prog
    
    p_b = np.array([(R + (r + offset) * np.cos(v_b)) * np.cos(u_b),
                    (R + (r + offset) * np.cos(v_b)) * np.sin(u_b),
                    (r + offset) * np.sin(v_b)])
    
    v_vec_b = np.array([-np.sin(v_b)*np.cos(u_b), -np.sin(v_b)*np.sin(u_b), np.cos(v_b)])
    v_vec_b = v_vec_b / np.linalg.norm(v_vec_b) # Normalize
    w_vec_b = np.array([-np.sin(u_b), np.cos(u_b), 0])
    
    tip_b = p_b + v_vec_b * 0.8
    base_b = tip_b - v_vec_b * 0.3
    left_b = base_b + w_vec_b * 0.2
    right_b = base_b - w_vec_b * 0.2
    
    # BLUE ARROW: Explicit bounds logic (-0.1*pi to 1.1*pi)
    # Modulo ensures the logic safely wraps around the 0 / 2pi boundary
    v_mod = v_b % (2 * np.pi)
    # -0.1*pi wrapped to modulo 2pi is 1.9*pi.
    is_blue_front = (v_mod <= 1.1 * np.pi) or (v_mod >= 1.9 * np.pi)
    blue_layer = 100 if is_blue_front else 0
    
    l3, = ax.plot([p_b[0], tip_b[0]], [p_b[1], tip_b[1]], [p_b[2], tip_b[2]], color="#284D78", lw=4, zorder=blue_layer)
    l4, = ax.plot([left_b[0], tip_b[0], right_b[0]], [left_b[1], tip_b[1], right_b[1]], [left_b[2], tip_b[2], right_b[2]], color="#284D78", lw=4, zorder=blue_layer)
    arrow_artists.extend([l3, l4])
    
# 4. Render and Save
parser = argparse.ArgumentParser(description='Generate toric loop animation GIF')
parser.add_argument('--dpi', type=int, default=100, help='DPI for the saved GIF (default: 100)')
args = parser.parse_args()

print(f"Rendering GIF at {args.dpi} DPI... This may take ~10 seconds.")
ani = FuncAnimation(fig, update, frames=np.arange(0, total_frames), interval=50)
ani.save('toric_loops.gif', writer=PillowWriter(fps=20), dpi=args.dpi)
print("Done! Check out 'toric_loops.gif'.")