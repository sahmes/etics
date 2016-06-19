namespace etics {
    namespace multicenter {
        void Init(int Nmax, int k3gs_new, int k3bs_new, int k4gs_new, int k4bs_new);
        void CalculateGravity(Particle *P, int N, Real *Potential, vec3 *F);
    }
}