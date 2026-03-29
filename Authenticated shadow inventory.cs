public class AuthenticatedShadowInventory
{
    private readonly ShadowInventory _shadow;
    private readonly InventoryArtifact _artifact;

    public AuthenticatedShadowInventory(ShadowInventory shadow, InventoryArtifact artifact)
    {
        _shadow = shadow;
        _artifact = artifact;
    }

    private bool CanStore(IItem item)
    {
        // Example: item has tags, artifact has allowed tags
        // If no tags, treat as restricted
        if (item.Tags == null || item.Tags.Count == 0)
            return false;

        return item.Tags.Any(tag => _artifact.AllowedItemTags.Contains(tag));
    }

    public bool AddItem(IItem item, int quantity)
    {
        if (_artifact.IsReadOnly)
            return false;

        if (!CanStore(item))
            return false;

        // Optional: enforce MaxSlots by counting distinct item ids
        // before delegating to shadow
        return _shadow.AddItem(item, quantity);
    }

    public bool RemoveItem(string itemId, int quantity)
    {
        if (_artifact.IsReadOnly)
            return false;

        return _shadow.RemoveItem(itemId, quantity);
    }

    public int GetQuantity(string itemId) => _shadow.GetQuantity(itemId);
}
