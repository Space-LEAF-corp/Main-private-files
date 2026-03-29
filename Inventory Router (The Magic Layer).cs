public class InventoryRouter : IInventory
{
    private readonly GameInventory gameInventory;
    private readonly ShadowInventory shadowInventory;

    public InventoryRouter(GameInventory gameInv, ShadowInventory shadowInv)
    {
        gameInventory = gameInv;
        shadowInventory = shadowInv;
    }

    public bool AddItem(IItem item, int quantity)
    {
        // Try game inventory first
        if (gameInventory.AddItem(item, quantity))
            return true;

        // Overflow goes to shadow inventory
        return shadowInventory.AddItem(item, quantity);
    }

    public bool RemoveItem(string itemId, int quantity)
    {
        // Try removing from game inventory first
        if (gameInventory.RemoveItem(itemId, quantity))
            return true;

        // If not enough, pull from shadow inventory
        return shadowInventory.RemoveItem(itemId, quantity);
    }

    public int GetQuantity(string itemId)
    {
        return gameInventory.GetQuantity(itemId)
             + shadowInventory.GetQuantity(itemId);
    }
}
