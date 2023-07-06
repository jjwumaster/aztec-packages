import { AztecNode, HttpNode } from '@aztec/aztec-node';
import { ContractDeployer, createAztecRPCServer, TxStatus, TxHash, Point, AztecAddress, Fr } from '@aztec/aztec.js';
import { DebugLogger } from '@aztec/foundation/log';
import { ContractAbi } from '@aztec/foundation/abi';
import { sleep } from '@aztec/foundation/sleep';
import { ZkTokenContractAbi } from '@aztec/noir-contracts/examples';
import { randomBytes } from 'crypto';

/**
 * Helper function for creating an instance of the aztec rpc server.
 * @param numberOfAccounts - The initial number of accounts to be added.
 * @param aztecNode - The instance of aztec node to be used by the rpc server.
 * @returns The instance of an aztec rpc server.
 */
export async function createAztecRpc(numberOfAccounts = 1, aztecNode: AztecNode) {
  const arc = await createAztecRPCServer(aztecNode);

  for (let i = 0; i < numberOfAccounts; ++i) {
    // Since we will be just deploying contracts, we don't really need an account
    // contract deployed, just a keypair. So we just register a random address as
    // the account contract implementation. This hack hints at the fact that we
    // need to rethink the APIs of the aztec-rpc-server and keystore post AA.
    const privKey = randomBytes(32);
    await arc.addAccount(privKey, AztecAddress.random(), Fr.random());
  }

  return arc;
}

const pointToPublicKey = (point: Point) => {
  const x = point.x.toBigInt();
  const y = point.y.toBigInt();
  return {
    x,
    y,
  };
};

/**
 * Helper function to pause until the provided aztec node instance is ready for use.
 * @param aztecNode - The instance of an aztec node that we need to wait for.
 * @param logger - An instance of a logging object.
 */
export async function waitUntilReady(aztecNode: AztecNode, logger: DebugLogger) {
  while (true) {
    try {
      const isReady = await aztecNode.isReady();
      if (!isReady) {
        throw new Error('Not ready');
      }
      break;
    } catch (err) {
      logger(`Aztec node not ready, will wait 10 seconds and check again...`);
      await sleep(10000);
    }
  }
}

/**
 * Function for carrying out the 'deployL2' command. An instance of the ZKToken Contract will be deployed
 * either once or repeatedly if specified.
 * @param rollupProviderUrl - The url of the rollup provider service the contracts should be deployed to.
 * @param intervalMs - The interval (ms) between repeated deployments. If 0, contract is only deployed once.
 * @param logger - An instance of a logging object.
 */
export async function deployL2Contract(rollupProviderUrl: string, intervalMs: number, logger: DebugLogger) {
  const node: AztecNode = new HttpNode(rollupProviderUrl);
  await waitUntilReady(node, logger);
  const aztecRpcServer = await createAztecRpc(2, node);
  const accounts = await aztecRpcServer.getAccounts();
  const outstandingTxs: TxHash[] = [];

  do {
    logger(`Deploying L2 contract...`);
    const initialBalance = 1_000_000_000n;
    const zkContract = ZkTokenContractAbi as ContractAbi;
    const owner = await aztecRpcServer.getAccountPublicKey(accounts[0]);
    const deployer = new ContractDeployer(zkContract, aztecRpcServer, owner);
    const d = deployer.deploy(initialBalance, pointToPublicKey(owner));
    await d.create();
    const tx = d.send();
    const txHash = await tx.getTxHash();
    outstandingTxs.push(txHash);
    logger('L2 contract deployed');
    await sleep(intervalMs > 0 ? intervalMs : 1);
    for (let i = 0; i < outstandingTxs.length; i++) {
      const hash = outstandingTxs[i];
      const receipt = await aztecRpcServer.getTxReceipt(hash);
      if (receipt.status == TxStatus.MINED) {
        logger(`Tx ${hash.toString()} settled, contract address ${receipt.contractAddress}`);
        outstandingTxs.splice(i, 1);
      }
    }
  } while (intervalMs != 0);
}
