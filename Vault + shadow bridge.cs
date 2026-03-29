public class WarehouseVault
{
    public string VaultId { get; set; }
    public Dictionary<VaultTier, VaultSection> Sections { get; set; } = new();
}

public class VaultBridge
{
    private readonly ShadowInventory _shadow;
    private readonly WarehouseVault _vault;
    private readonly string _playerId;

    public VaultBridge(ShadowInventory shadow, WarehouseVault vault, string playerId)
    {
        _shadow = shadow;
        _vault = vault;
        _playerId = playerId;
    }

    public bool SendToVault(string itemId, int quantity, VaultTier tier)
    {
        if (_shadow.GetQuantity(itemId) < quantity)
            return false;

        _shadow.RemoveItem(itemId, quantity);

        var section = _vault.Sections.GetValueOrDefault(tier)
                      ?? (_vault.Sections[tier] = new VaultSection { Tier = tier });

        if (!section.Slots.TryGetValue(itemId, out var slot))
        {
            slot = new VaultSlot
            {
                ItemId = itemId,
                Quantity = 0,
                StoredAt = DateTime.UtcNow,
                StoredByPlayerId = _playerId,
                Provenance = "ShadowInventory"
            };
            section.Slots[itemId] = slot;
        }

        slot.Quantity += quantity;
        return true;
    }

    public bool RetrieveFromVault(string itemId, int quantity, VaultTier tier, IItem itemDef)
    {
        var section = _vault.Sections.GetValueOrDefault(tier);
        if (section == null || !section.Slots.TryGetValue(itemId, out var slot))
            return false;

        if (slot.Quantity < quantity)
            return false;

        if (!_shadow.AddItem(itemDef, quantity))
            return false;

        slot.Quantity -= quantity;
        return true;
    }
}
