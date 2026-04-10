// katana-duel-server.ts
// Minimal, event-driven duel coordinator. No external IP; all generic “katana” semantics.

import { WebSocketServer } from "ws";

type Player = { id: string; socket: any; activeKatana: boolean; inDuel: boolean };
type Duel = { id: string; p1: Player; p2: Player; state: "pending" | "active" | "ended" };

const wss = new WebSocketServer({ port: 8081 });
const players = new Map<string, Player>();
const queue: Player[] = [];
const duels = new Map<string, Duel>();

function broadcastDuelState(duel: Duel) {
  const payload = JSON.stringify({ type: "DUEL_STATE", duelId: duel.id, state: duel.state });
  duel.p1.socket.send(payload);
  duel.p2.socket.send(payload);
}

function tryMatch() {
  while (queue.length >= 2) {
    const p1 = queue.shift()!;
    const p2 = queue.shift()!;
    const duelId = `duel_${Date.now()}_${Math.random().toString(36).slice(2)}`;
    const duel = { id: duelId, p1, p2, state: "active" as const };
    p1.inDuel = true; p2.inDuel = true;
    duels.set(duelId, duel);
    broadcastDuelState(duel);
  }
}

wss.on("connection", (socket) => {
  const id = `p_${Date.now()}_${Math.random().toString(36).slice(2)}`;
  const player: Player = { id, socket, activeKatana: false, inDuel: false };
  players.set(id, player);

  socket.send(JSON.stringify({ type: "WELCOME", playerId: id }));

  socket.on("message", (raw: string) => {
    const msg = JSON.parse(raw);
    switch (msg.type) {
      case "SET_ACTIVE_KATANA":
        player.activeKatana = !!msg.active;
        // Only enqueue players who opt in and are not already dueling
        if (player.activeKatana && !player.inDuel && !queue.includes(player)) queue.push(player);
        tryMatch();
        break;

      case "DUEL_INPUT":
        // Forward minimal input (e.g., swing events) to opponent during active duel
        const duel = [...duels.values()].find(d => d.p1.id === player.id || d.p2.id === player.id);
        if (!duel || duel.state !== "active") return;
        const opponent = duel.p1.id === player.id ? duel.p2 : duel.p1;
        opponent.socket.send(JSON.stringify({ type: "OPPONENT_INPUT", data: msg.data }));
        break;

      case "END_DUEL":
        {
          const d = duels.get(msg.duelId);
          if (!d) return;
          d.state = "ended";
          d.p1.inDuel = false; d.p2.inDuel = false;
          broadcastDuelState(d);
          duels.delete(msg.duelId);
          // Re-queue players if they still want katana active
          [d.p1, d.p2].forEach(p => { if (p.activeKatana) queue.push(p); });
          tryMatch();
        }
        break;
    }
  });

  socket.on("close", () => {
    players.delete(id);
    const idx = queue.indexOf(player);
    if (idx >= 0) queue.splice(idx, 1);
    // Clean up duels
    const d = [...duels.values()].find(x => x.p1.id === id || x.p2.id === id);
    if (d) {
      d.state = "ended";
      broadcastDuelState(d);
      duels.delete(d.id);
    }
  });
});
// DualWieldController.cs
// Attach to player avatar. Supports sword+sword or sword+shield.

using UnityEngine;

public enum WeaponType { None, Sword, Shield }

[System.Serializable]
public class WeaponSlot {
    public string name;
    public WeaponType type;
    public Transform weaponTransform;
    public Collider weaponCollider;
}

public class DualWieldController : MonoBehaviour
{
    [Header("Weapon Slots")]
    public WeaponSlot rightHand;
    public WeaponSlot leftHand;

    [Header("Combat Settings")]
    public AudioSource audioSource;
    public AudioClip[] metalClashClips;
    public ParticleSystem sparkVFX;
    public LayerMask bladeCollisionMask;

    [Header("Shield Settings")]
    public float blockAngle = 120f; // degrees of effective block arc
    public float blockStaminaCost = 10f;

    private bool isBlocking;

    void Update()
    {
        HandleInput();
    }

    void HandleInput()
    {
        // Example input scheme: Fire1 = right hand, Fire2 = left hand, Block = shield
        if (Input.GetButtonDown("Fire1") && rightHand.type == WeaponType.Sword)
            Swing(rightHand);

        if (Input.GetButtonDown("Fire2") && leftHand.type == WeaponType.Sword)
            Swing(leftHand);

        if (leftHand.type == WeaponType.Shield)
        {
            isBlocking = Input.GetButton("Block");
            if (isBlocking) RaiseShield();
        }
    }

    void Swing(WeaponSlot slot)
    {
        // Trigger swing animation for that hand
        Animator anim = GetComponent<Animator>();
        if (anim) anim.SetTrigger(slot.name + "_Swing");

        // Trace collisions for clash
        TraceWeaponForClash(slot);
    }

    void TraceWeaponForClash(WeaponSlot slot)
    {
        if (!slot.weaponTransform || !slot.weaponCollider) return;
        Collider[] hits = Physics.OverlapSphere(slot.weaponTransform.position, 0.1f, bladeCollisionMask);
        foreach (var h in hits)
        {
            if (h.transform.root != this.transform.root)
            {
                PlayMetalClashSFX();
                EmitSparksAt(slot.weaponTransform.position, h.transform);
            }
        }
    }

    void RaiseShield()
    {
        // Simple block mechanic: check incoming attack angle
        // In practice, integrate with opponent swing direction
        // Here we just set a blocking state
        Debug.Log("Shield raised: blocking active");
    }

    void PlayMetalClashSFX()
    {
        if (!audioSource || metalClashClips.Length == 0) return;
        var clip = metalClashClips[Random.Range(0, metalClashClips.Length)];
        audioSource.PlayOneShot(clip);
    }

    void EmitSparksAt(Vector3 point, Transform target)
    {
        if (!sparkVFX) return;
        sparkVFX.transform.position = point;
        sparkVFX.transform.LookAt(target.position, Vector3.up);
        sparkVFX.Play();
    }
}
