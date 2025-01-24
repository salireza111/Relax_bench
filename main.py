import pygame
import sys
import math
import threading
import time
import copy
import platform
import subprocess
import psutil
import numpy as np
import pyrosetta
import cpuinfo
import os
from pygame.locals import *
from pyrosetta import Pose, get_fa_scorefxn
from pyrosetta.rosetta.core.scoring import CA_rmsd, hbonds
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.scoring.sasa import SasaCalc

WIDTH, HEIGHT = 1280, 720
FONT_SIZE = 18
COLORS = {
    'background': (30, 30, 30),
    'button': (50, 200, 50),
    'hover': (70, 220, 70),
    'text': (255, 255, 255),
    'backbone': (200, 25, 25),      
    'other_atoms': (100, 100, 255),   
    'bond': (10, 10, 10),             
    'panel_bg': (40, 40, 40),
    'panel_border': (200, 200, 200),
    'button_border': (180, 180, 180),
    'gradient_start': (50, 50, 50),
    'gradient_end': (100, 100, 100),
    'button_highlight': (100, 250, 100),
}

def get_cpu_info():
    try:
        info = cpuinfo.get_cpu_info()
        return (info['brand_raw'])
    except:
        return platform.processor()

def calculate_hbonds(pose):
    hbonds_set = hbonds.HBondSet()
    pose.update_residue_neighbors()
    hbonds.fill_hbond_set(pose, False, hbonds_set)
    return hbonds_set.nhbonds()

def parse_pose(pose):
    atoms = []
    backbone_bonds = []
    residue_map = {}
    
    for resi in range(1, pose.total_residue() + 1):
        residue = pose.residue(resi)
        res_key = (pose.pdb_info().chain(resi), str(resi)) if pose.pdb_info() else ('A', str(resi))
        residue_map[res_key] = {'N': None, 'CA': None, 'C': None}
        
        for atom in range(1, residue.natoms() + 1):
            atom_name = residue.atom_name(atom).strip()
            xyz = residue.xyz(atom)
            idx = len(atoms)
            atoms.append({
                'x': xyz.x, 'y': xyz.y, 'z': xyz.z,
                'residue': res_key, 'atom_name': atom_name
            })
            
            if atom_name in ['N', 'CA', 'C']:
                residue_map[res_key][atom_name] = idx

    sorted_res = sorted(residue_map.keys(), key=lambda x: (x[0], int(x[1])))
    for i, res in enumerate(sorted_res):
        if 'N' in residue_map[res] and 'CA' in residue_map[res]:
            backbone_bonds.append((residue_map[res]['N'], residue_map[res]['CA']))
        if 'CA' in residue_map[res] and 'C' in residue_map[res]:
            backbone_bonds.append((residue_map[res]['CA'], residue_map[res]['C']))
        if i < len(sorted_res)-1:
            next_res = sorted_res[i+1]
            if 'C' in residue_map[res] and 'N' in residue_map[next_res]:
                backbone_bonds.append((residue_map[res]['C'], residue_map[next_res]['N']))
    
    return atoms, backbone_bonds

def rotate_point(x, y, z, angle_x, angle_y):
    rad_y = math.radians(angle_y)
    x_rot = x * math.cos(rad_y) + z * math.sin(rad_y)
    z_rot = z * math.cos(rad_y) - x * math.sin(rad_y)
    
    rad_x = math.radians(angle_x)
    y_rot = y * math.cos(rad_x) - z_rot * math.sin(rad_x)
    z_final = y * math.sin(rad_x) + z_rot * math.cos(rad_x)
    
    return x_rot, y_rot, z_final

def project_point(x, y, z, width, height, fov, dist):
    z_clip = z + dist
    if z_clip <= 0: z_clip = 0.001
    factor = fov / z_clip
    return (
        int(x * factor + width/2),
        int(-y * factor + height/2)
    )

class RelaxationController:
    def __init__(self, pose, atoms, lock):
        self.pose = pose
        self.initial_pose = pose.clone()
        self.atoms = atoms
        self.lock = lock
        self.running = False
        self.completed = False

        self.metrics = {
            'time': {'total': 0, 'per_iter': []},
            'memory': {'max': 0, 'usage': []},
            'hbonds': {'initial': 0, 'final': 0},
            'sasa': {'initial': 0, 'final': 0},
            'energy': {'initial': 0, 'final': 0},
            'rmsd': 0,
            'residues': pose.total_residue(),
            'atoms': len(atoms)
        }

        self._collect_initial_metrics()

    def _collect_initial_metrics(self):
        self.metrics['energy']['initial'] = get_fa_scorefxn().score(self.initial_pose)
        self.metrics['hbonds']['initial'] = calculate_hbonds(self.initial_pose)
        sasa_calc = SasaCalc()
        sasa_calc.calculate(self.initial_pose)
        self.metrics['sasa']['initial'] = sasa_calc.get_total_sasa()

    def check_convergence(self, energy_initial, energy_current, energy_threshold=0.1):
        return abs(energy_current - energy_initial) < energy_threshold

    def start_relaxation(self):
        if not self.running:
            self.running = True
            self.thread = threading.Thread(target=self.run_relaxation)  
            self.thread.start()

    def run_relaxation(self): 
        start_time = time.time()
        start_mem = psutil.Process().memory_info().rss

        scorefxn = get_fa_scorefxn()
        mm = MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)

        min_mover = MinMover()
        min_mover.movemap(mm)
        min_mover.score_function(scorefxn)
        min_mover.max_iter(8)

        energy_initial = scorefxn.score(self.pose)

        for loop in range(120):  
            if not self.running:
                break

            iter_start = time.time()
            min_mover.apply(self.pose)

            energy_current = scorefxn.score(self.pose)

            if self.check_convergence(energy_initial, energy_current, energy_threshold=0.1):
                print(f"Relaxation stopped at loop {loop} due to minimal energy change.")
                break  

            energy_initial = energy_current  

            with self.lock:
                for atom in self.atoms:
                    res_num = int(atom['residue'][1])
                    try:
                        residue = self.pose.residue(res_num)
                        atom_idx = residue.atom_index(atom['atom_name'])
                        xyz = residue.xyz(atom_idx)
                        atom.update(x=xyz.x, y=xyz.y, z=xyz.z)
                    except:
                        continue

            self.metrics['time']['per_iter'].append(time.time() - iter_start)
            self.metrics['memory']['usage'].append(
                psutil.Process().memory_info().rss - start_mem
            )
            time.sleep(0.01)

        self.metrics['time']['total'] = time.time() - start_time
        self.metrics['memory']['max'] = max(self.metrics['memory']['usage'])
        self.metrics['energy']['final'] = scorefxn.score(self.pose)
        self.metrics['hbonds']['final'] = calculate_hbonds(self.pose)
        self.metrics['rmsd'] = CA_rmsd(self.initial_pose, self.pose)

        sasa_calc = SasaCalc()
        sasa_calc.calculate(self.pose)
        self.metrics['sasa']['final'] = sasa_calc.get_total_sasa()

        self.completed = True
        self.running = False


class PDBViewer:
    def __init__(self, pdb_file):
        pygame.init()
        self.screen = pygame.display.set_mode((WIDTH, HEIGHT))
        pygame.display.set_caption("Protein Relaxation Benchmark")
        self.clock = pygame.time.Clock()
        self.font = pygame.font.SysFont('Arial', FONT_SIZE)
        
        pyrosetta.init()
        self.pose = Pose()
        pyrosetta.pose_from_file(self.pose, pdb_file)
        
        self.data_lock = threading.Lock()
        self.atoms, self.bonds = parse_pose(self.pose)
        self.controller = RelaxationController(self.pose, self.atoms, self.data_lock)
        
        self.view = {
            'angle_x': 0,
            'angle_y': 0,
            'distance': 300,
            'fov': 1000,
            'dragging': False
        }
        
        self.cpu_info = get_cpu_info()

    def run(self):
        while True:
            self.handle_events()
            self.draw()
            self.clock.tick(60)

    def handle_events(self):
        mouse_pos = pygame.mouse.get_pos()
        for event in pygame.event.get():
            if event.type == QUIT:
                pygame.quit()
                sys.exit()
                
            if event.type == MOUSEBUTTONDOWN:
                if event.button == 1:
                    self.view['dragging'] = True
                    pygame.mouse.get_rel()
                    
                    if not self.controller.running and not self.controller.completed:
                        btn_rect = pygame.Rect(WIDTH-220, HEIGHT-80, 200, 40)
                        if btn_rect.collidepoint(mouse_pos):
                            self.controller.start_relaxation()
                            
            elif event.type == MOUSEBUTTONUP:
                if event.button == 1:
                    self.view['dragging'] = False
                    
            elif event.type == MOUSEMOTION and self.view['dragging']:
                dx, dy = pygame.mouse.get_rel()
                self.view['angle_y'] += dx * 0.3
                self.view['angle_x'] += dy * 0.3
                
            elif event.type == MOUSEWHEEL:
                self.view['distance'] = max(100, self.view['distance'] - event.y * 20)
    

    def draw(self):
        self.screen.fill(COLORS['background'])
        self.draw_gradient()
        
        with self.data_lock:
            atoms_copy = copy.deepcopy(self.atoms)
        
        projected = []
        for idx, atom in enumerate(atoms_copy):
            x, y, z = rotate_point(
                atom['x'], atom['y'], atom['z'],
                self.view['angle_x'],
                self.view['angle_y']
            )
            px, py = project_point(
                x, y, z, WIDTH, HEIGHT,
                self.view['fov'], self.view['distance']
            )
            depth = z + self.view['distance']  
            projected.append((px, py, depth, idx))  

        for i, j in self.bonds:
            if i < len(projected) and j < len(projected):
                px1, py1, _, _ = projected[i]
                px2, py2, _, _ = projected[j]
                pygame.draw.line(self.screen, COLORS['bond'], (px1, py1), (px2, py2), 1)

        projected_sorted = sorted(projected, key=lambda x: -x[2])

        for px, py, _, idx in projected_sorted:
            atom = atoms_copy[idx]
            color = COLORS['backbone'] if atom['atom_name'] in ['N', 'CA', 'C'] \
                   else COLORS['other_atoms']
            radius = 4 if atom['atom_name'] in ['N', 'CA', 'C'] else 2.5
            pygame.draw.circle(self.screen, color, (px, py), int(radius))

        self.draw_ui_elements()
        pygame.display.flip()


    def draw_gradient(self):
        for i in range(HEIGHT):
            color = tuple(
                int(COLORS['gradient_start'][x] + (COLORS['gradient_end'][x] - COLORS['gradient_start'][x]) * (i / HEIGHT)) 
                for x in range(3)
            )
            pygame.draw.line(self.screen, color, (0, i), (WIDTH, i))

    def draw_ui_elements(self):
        cpu_text = self.font.render(f"CPU: {self.cpu_info}", True, COLORS['text'])
        self.screen.blit(cpu_text, (20, 20))
        
        if not self.controller.completed:
            if not self.controller.running:
                btn_rect = pygame.Rect(WIDTH-220, HEIGHT-80, 200, 40)
                self.draw_button(WIDTH-220, HEIGHT-80, 200, 40, "Start Benchmark")
            else:
                self.draw_text(WIDTH-220, HEIGHT-80, "Running...", COLORS['text'])

        if self.controller.completed:
            self.draw_results_panel()

    def draw_button(self, x, y, w, h, text):
        btn_rect = pygame.Rect(x, y, w, h)
        color = COLORS['hover'] if btn_rect.collidepoint(pygame.mouse.get_pos()) \
               else COLORS['button']
        pygame.draw.rect(self.screen, color, btn_rect, border_radius=10)
        pygame.draw.rect(self.screen, COLORS['button_border'], btn_rect, 3, border_radius=10)
        text_surf = self.font.render(text, True, (255,255,255))
        text_rect = text_surf.get_rect(center=btn_rect.center)
        self.screen.blit(text_surf, text_rect)

    def draw_text(self, x, y, text, color):
        text_surf = self.font.render(text, True, color)
        self.screen.blit(text_surf, (x, y))

    def draw_results_panel(self):
        panel_width = 400
        panel = pygame.Surface((panel_width, 200))
        panel.fill(COLORS['panel_bg'])
        pygame.draw.rect(panel, COLORS['panel_border'], panel.get_rect(), 2, border_radius=8)
        
        delta_energy = self.controller.metrics['energy']['final'] - self.controller.metrics['energy']['initial']
        total_time = self.controller.metrics['time']['total']

        score = (22.5 - (delta_energy / total_time)) * 1000

        metrics = [
            f"Residues: {self.controller.metrics['residues']}",
            f"Rlax_Bench Score: {score:.2f}",
            f"Time: {self.controller.metrics['time']['total']:.2f}s",
            f"ΔEnergy: {delta_energy:.2f}",
            f"ΔHBonds: {self.controller.metrics['hbonds']['final'] - self.controller.metrics['hbonds']['initial']}",
            f"ΔSASA: {self.controller.metrics['sasa']['final'] - self.controller.metrics['sasa']['initial']:.1f}",
            f"Cα RMSD: {self.controller.metrics['rmsd']:.3f}Å",
            f"Score: {score:.2f}"
        ]

        
        y = 20
        for line in metrics:
            text = self.font.render(line, True, COLORS['text'])
            panel.blit(text, (20, y))
            y += 25
        
        self.screen.blit(panel, (WIDTH - panel_width - 20, 20))

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))
    viewer = PDBViewer(script_dir + "/RNR.pdb")  
    viewer.run()

#132 -1280
#237 1519